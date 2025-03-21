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

__all__ = ["SingleFrameScoreMeasurementTask"]

import numpy as np
# import scipy.signal

import lsst.afw.image
import lsst.geom
import lsst.meas.base


class SingleFrameScoreMeasurementConfig(lsst.meas.base.BaseMeasurementConfig):
    """Config class for single frame measurement driver task.
    """

    plugins = lsst.meas.base.SingleFramePlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to be run and their configuration"
    )
    undeblended = lsst.meas.base.SingleFramePlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to run on undeblended image"
    )

import lsst.afw.display
display = lsst.afw.display.Display()
display.frame = 3


class SingleFrameScoreMeasurementTask(lsst.meas.base.SingleFrameMeasurementTask):
    """Run measurement plugins on a likelihood/score image, like that produced
    by `~lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask`.
    """
    # TODO: need to create an appropriate "cannot deconvolve" flag
    # TODO: will have to override run and runPlugins, since I need to
    # pass other arguments.

    # For FFT deconvolve, the kernel and cutout must be the same size.
    # Can we get a gaussian parameterization of the PSF instead of using
    # a PSF image? If so, then calculate the image of the PSF within the
    # cutout box, and use that.
    # use np.fft for simplicity

    def callMeasure(self, record, exposure, kernel, *args, **kwargs):
        beginOrder = kwargs.pop("beginOrder", None)
        endOrder = kwargs.pop("endOrder", None)

        self.doMeasurement(self.plugins["base_SdssCentroid"], record, exposure=exposure,
                           *args, **kwargs)

        # TODO: either we have to increase the size of the kernel image box, or
        # the size of the cutout, so that they match. Always use the largest size.
        box = record.getFootprint().getBBox()
        # padded = lsst.geom.Box2I(box)
        # padded.grow(box.height//2)
        dim = box.getDimensions()
        dim.x += 1
        dim.y += 1
        box = lsst.geom.Box2I.makeCenteredBox(box.getCenter(), dim)
        kernelBox = kernel.computeBBox(record.getCentroid())
        # if (box.width < kernelBox.width) or (box.height < kernelBox.height):
        #     raise RuntimeError("Can't handle this yet!")

        # center = lsst.geom.Point2D(record.getCentroid())
        # center.shift(lsst.geom.Extent2D(1, 1))
        # localKernel = kernel.getLocalKernel(record.getCentroid())
        kernelImage = lsst.afw.image.ImageD(box)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        tempBox = lsst.geom.Box2I.makeCenteredBox(box.getCenter(), kernelBox.getDimensions())
        # localKernel.computeImage(kernelImage[tempBox], doNormalize=True)
        kernelImage[tempBox] = kernel.computeKernelImage(record.getCentroid())
        # kernelImage[box.getCenter()] = 1
        deconvolved = lsst.afw.image.ExposureD(box)
        shape = deconvolved.image.array.shape
        cutout = lsst.afw.image.ImageF(box)
        cutout[box].array = exposure.image[box].array
        # TODO: restore rfft2 once we've debugged things (will make slicing harder)
        fftcutout = np.fft.fft2(exposure.image[box].array, shape)
        fftkernel = np.fft.fft2(np.roll(kernelImage.array, (0, 0)), shape)
        # TODO: consider a blackman-harris windowing function on either the
        # quotient or fftcutout.
        quotient = np.zeros_like(fftkernel)
        # TODO: How do we choose the best region to remove?
        slice = np.s_[12:17, 12:17]
        quotient[slice] = np.fft.fftshift(fftcutout / np.abs(fftkernel))[slice]
        deconvolved.image.array = np.fft.ifft2(np.fft.fftshift(quotient),
                                               deconvolved.image.array.shape).real
        temp = lsst.afw.image.ImageD(box)
        display.frame += 1
        display.image(kernelImage, title="kernelImage")
        display.frame += 1
        display.image(deconvolved, title="deconvolved")
        display.frame += 1
        display.image(cutout, title="exposure")
        display.frame += 1
        temp.array = np.fft.fftshift(np.array(np.abs(fftcutout), dtype=np.float64))
        display.image(temp, title="fftcutout")
        display.frame += 1
        temp.array = np.fft.fftshift(np.abs(fftkernel))
        display.image(temp, title="fftkernel")
        display.frame += 1
        temp.array = np.abs(quotient)
        display.image(temp, title="fftcutout/fftkernel")

        # display.frame = 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (0, 0))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="no shift")
        # display.frame += 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (1, 0))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="x=1")
        # display.frame += 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (-1, 0))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="x=-1")
        # display.frame += 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (0, 1))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="y=1")
        # display.frame += 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (0, -1))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="y=-1")
        # display.frame += 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (1, 1))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="both=1")
        # display.frame += 1
        # fftkernel = np.fft.fft2(np.roll(kernelImage.array, (-1, -1))))
        # deconvolved.image.array = np.fft.ifft2(fftcutout/fftkernel),
        #                                        deconvolved.image.array.shape).real
        # display.image(deconvolved, title="both=-1")

        # # display.frame = 1
        # # display.image(fftcutout, title="fftcutout")
        # # display.frame = 1
        # # display.image(fftkernel, title="fftkernel")
        # # display.frame = 1
        # # display.image(fftcutout/fftkernel, title="fftcutout/fftkernel")

        import os; print(os.getpid()); import ipdb; ipdb.set_trace();

        popped = self.plugins.pop("base_SdssCentroid")
        for plugin in self.plugins.iter():
            if beginOrder is not None and plugin.getExecutionOrder() < beginOrder:
                continue
            if endOrder is not None and plugin.getExecutionOrder() >= endOrder:
                break
            self.doMeasurement(plugin, record, exposure=deconvolved,
                               *args, **kwargs)

        self.plugins["base_SdssCentroid"] = popped
        return
