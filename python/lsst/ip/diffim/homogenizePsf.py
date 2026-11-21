# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["HomogenizePsfConfig", "HomogenizePsfTask"]

import numpy as np

import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod

from .utils import _interpolateImage, sampleCoaddPsfMoments


class HomogenizePsfConfig(pexConfig.Config):
    """Configuration for `HomogenizePsfTask`."""

    quantile = pexConfig.RangeField(
        dtype=float,
        default=0.8,
        min=0.0,
        max=1.0,
        inclusiveMax=True,
        doc="Fraction of the template area whose PSF is at most as wide as the"
        " common width. At most the remaining fraction is dropped; 1 drops"
        " nothing and broadens everything to the widest PSF.",
    )
    tolerance = pexConfig.RangeField(
        dtype=float,
        default=0.03,
        min=0.0,
        doc="Fractional amount by which the PSF width of a part of the"
        " template may exceed the common width and still be kept, unchanged."
        " PSF matching absorbs mismatches this small.",
    )
    gridStep = pexConfig.RangeField(
        dtype=int,
        default=24,
        min=1,
        doc="Size in pixels of the square blocks of the template that are each"
        " broadened to the common width as a unit. Where the PSF steps inside"
        " a block, part of the block gets the wrong broadening, so smaller"
        " blocks leave less residual mismatch but take longer to measure.",
    )
    broadeningStep = pexConfig.RangeField(
        dtype=float,
        default=0.5,
        min=0.0,
        inclusiveMin=False,
        doc="Sigma, in pixels, of the narrowest non-zero Gaussian kernel used to"
        " convolve the template. Level n of the ladder of kernels has sigma"
        " broadeningStep*sqrt(n), so that their sigma^2 are evenly spaced. Each"
        " block receives the convolved templates of the two kernels that bracket"
        " the Gaussian matching kernel that would give a uniform PSF size,"
        " linearly interpolated in sigma^2.",
    )
    badMaskPlanes = pexConfig.ListField(
        dtype=str,
        default=("NO_DATA", "BAD", "SAT", "EDGE"),
        doc="Mask planes whose pixels are filled before the template is"
        " convolved, so that their values do not spread into the pixels"
        " around them. They keep their original values afterwards.",
    )


class HomogenizePsfTask(pipeBase.Task):
    """Broaden the PSF of a template to a common width, and drop the parts of
    the template where the PSF is much wider.

    A template assembled from coadds has a PSF that steps wherever the set of
    input images changes, which PSF matching cannot follow. A homogenized
    template has the same PSF everywhere it has data.
    """

    ConfigClass = HomogenizePsfConfig
    _DefaultName = "homogenizePsf"

    @timeMethod
    def run(self, template, psf, kernelSigmas=5):
        """Broaden the template PSF to a common width, and drop the parts of
        the template whose PSF is wider.

        The width of the CoaddPsf is measured on a grid of
        ``config.gridStep`` pixel blocks. The common width is the block width
        at ``config.quantile`` of the blocks. Each narrower block is convolved
        with the circular Gaussian whose sigma^2 raises the PSF's second
        moment to the common one. Blocks wider than the common width by more
        than ``config.tolerance`` are set to ``NO_DATA``, and those within the
        tolerance are kept unchanged. The template is convolved with a ladder
        of Gaussian kernels whose sigma^2 are spaced by the square of
        ``config.broadeningStep``, and each pixel is linearly interpolated
        between the two that bracket the broadening (added sigma^2) its block
        needs, which gives that block the needed second moment. Ellipticity
        and the spatial variation of the PSF inside a block are ignored.

        Only the image is convolved. Pixels in ``config.badMaskPlanes`` are
        filled before convolving, so that they do not spread into their
        neighbors, and get their original values back afterwards. The
        variance plane is left unchanged: convolving with a normalized kernel
        lowers the noise variance of each pixel but correlates neighboring
        pixels, and leaves the noise on scales larger than the kernel, which
        is what detection on the PSF-smoothed difference image measures, as
        it was.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            Template to homogenize in place. The image of the dropped blocks is
            set to 0, their variance plane to 0 and their ``NO_DATA`` mask
            plane, like the pixels that no coadd covers.
        psf : `lsst.meas.algorithms.CoaddPsf`
            PSF of the template, built from its inputs. Its coadd WCS is the
            template WCS.
        kernelSigmas : `float`, optional
            Half-width, in sigma, of the Gaussian kernels the template is
            convolved with.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Results as a struct with attributes:

            ``psfSigma``
                Gaussian sigma in pixels of the homogenized template PSF, or
                `None` if the PSF cannot be evaluated anywhere and the
                template was left unchanged (`float` or `None`).
        """
        bbox = template.getBBox()
        step = self.config.gridStep
        nx = -(-bbox.width//step)
        ny = -(-bbox.height//step)
        xc = np.minimum(bbox.minX + step*(np.arange(nx) + 0.5), bbox.maxX)
        yc = np.minimum(bbox.minY + step*(np.arange(ny) + 0.5), bbox.maxY)
        xg, yg = np.meshgrid(xc, yc)
        coaddWcs = psf.getCoaddWcs()
        ra, dec = coaddWcs.pixelToSkyArray(xg.ravel(), yg.ravel())
        ixx, iyy, _ = sampleCoaddPsfMoments(psf, ra, dec, coaddWcs)
        # Squared width (sigma^2) of the circular Gaussian with the same trace.
        width2 = ((ixx + iyy)/2).reshape(ny, nx)
        good = np.isfinite(width2)
        if not good.any():
            self.log.warning("Not homogenizing the template PSF: it cannot be evaluated anywhere.")
            return pipeBase.Struct(psfSigma=None)
        # The width of an actual block, not one interpolated between blocks.
        target2 = float(np.quantile(width2[good], self.config.quantile, method="lower"))
        tooWide = good & (width2 > target2*(1 + self.config.tolerance)**2)
        need = np.where(good, np.clip(target2 - width2, 0, None), 0.0)

        # Expand the block values to pixels.
        def toPixels(blocks):
            return np.repeat(np.repeat(blocks, step, axis=0), step, axis=1)[:bbox.height, :bbox.width]
        need = toPixels(need)
        drop = toPixels(tooWide)

        dSigma2 = self.config.broadeningStep**2
        nLevels = int(np.ceil(need.max()/dSigma2)) + 1
        filled = template.maskedImage.clone()
        _interpolateImage(filled, self.config.badMaskPlanes)
        badPlanes = [plane for plane in self.config.badMaskPlanes
                     if plane in filled.mask.getMaskPlaneDict()]
        bad = (filled.mask.array & filled.mask.getPlaneBitMask(badPlanes)) != 0
        image = np.zeros_like(filled.image.array)
        control = afwMath.ConvolutionControl()
        control.setDoNormalize(True)
        control.setDoCopyEdge(True)
        for level in range(nLevels):
            # Linear interpolation weight of this level at each pixel.
            weight = np.clip(1 - np.abs(need/dSigma2 - level), 0, None)
            if not weight.any():
                continue
            if level == 0:
                convolved = filled.image
            else:
                sigma = np.sqrt(level*dSigma2)
                size = 2*int(np.ceil(kernelSigmas*sigma)) + 1
                gaussian = afwMath.GaussianFunction1D(sigma)
                kernel = afwMath.SeparableKernel(size, size, gaussian, gaussian)
                convolved = filled.image.clone()
                afwMath.convolve(convolved, filled.image, kernel, control)
            image += weight*convolved.array
            del convolved
        image[bad] = template.image.array[bad]
        # Dropped blocks look like pixels that no coadd covers.
        image[drop] = 0.0
        template.image.array[:] = image
        template.variance.array[drop] = 0.0
        template.mask.array[drop] |= template.mask.getPlaneBitMask("NO_DATA")

        sigmaTarget = float(np.sqrt(target2))
        droppedFraction = float(tooWide.sum()/good.sum())
        self.metadata["broadenMax"] = float(need.max())
        self.metadata["droppedFraction"] = droppedFraction
        self.log.info("Homogenized the template PSF to a Gaussian of sigma %.3f pixels, broadening by up to"
                      " %.3f pixels^2 and dropping %.1f%% of the template whose PSF is wider.",
                      sigmaTarget, need.max(), 100*droppedFraction)
        return pipeBase.Struct(psfSigma=sigmaTarget)
