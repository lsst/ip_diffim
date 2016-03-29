from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

from collections import namedtuple
import numpy as np

# LSST imports
from lsst.afw.geom import Point2D
from lsst.afw.image import (ImageF, MaskedImageF, PARENT)

__all__ = ("DipoleFitAlgorithm")


class DipoleFitAlgorithm(object):
    """
    Lightweight class containing methods for fitting dipoles in diffims, used by DipoleFitPlugin.
    This code is documented in DMTN-007.
    """
    # Below is a (somewhat incomplete) list of improvements
    # that would be worth investigating, given the time:

    # 1. Initial fast test whether a background gradient needs to be fit
    # 2. Initial fast estimate of background gradient(s) params -- perhaps using numpy.lstsq
    # 3. better estimate for staring flux when there's a strong gradient
    # 4. evaluate necessity for separate parameters for pos- and neg- images
    # 5. only fit background OUTSIDE footprint and dipole params INSIDE footprint?
    # 6. requires a new package `lmfit` -- investiate others? (astropy/scipy/iminuit?)
    # 7. account for PSFs that vary across the exposures

    # Only import lmfit if someone wants to use the new DipoleFitAlgorithm.
    import lmfit

    # This is just a private version number to sync with the ipython notebooks that I have been
    # using for algorithm development.
    __private_version__ = '0.0.3'

    # Create a namedtuple to hold all of the relevant output from the lmfit results
    resultsOutput = namedtuple('resultsOutput',
                               ['psfFitPosCentroidX', 'psfFitPosCentroidY',
                                'psfFitNegCentroidX', 'psfFitNegCentroidY', 'psfFitPosFlux',
                                'psfFitNegFlux', 'psfFitPosFluxSigma', 'psfFitNegFluxSigma',
                                'psfFitCentroidX', 'psfFitCentroidY', 'psfFitOrientation',
                                'psfFitSignaltoNoise', 'psfFitChi2', 'psfFitRedChi2'])

    def __init__(self, diffim, posImage=None, negImage=None):
        """ Algorithm to run dipole measurement on a diaSource

        Parameters
        ----------
        diffim : afw.image.Exposure
           Exposure on which the diaSources were detected
        posImage : afw.image.Exposure
           "Positive" exposure from which the template was subtracted
        negImage : afw.image.Exposure
           "Negative" exposure which was subtracted from the posImage
        """

        self.diffim = diffim
        self.posImage = posImage
        self.negImage = negImage

    @staticmethod
    def genBgGradientModel(in_x, b=None, x1=0., y1=0., xy=None, x2=0., y2=0.):
        """Generate gradient model (2-d array) with up to 2nd-order polynomial

        Parameters
        ----------
        in_x : 3-d numpy.array
           Provides the input x,y grid upon which to compute the gradient
        b : float
           Intercept (0th-order parameter) for gradient. If None, do nothing (for speed).
        x1 : float, optional
           X-slope (1st-order parameter) for gradient.
        y1 : float, optional
           Y-slope (1st-order parameter) for gradient.
        xy : float, optional
           X*Y coefficient (2nd-order parameter) for gradient.
           If None, do not compute 2nd order polynomial.
        x2 : float, optional
           X**2 coefficient (2nd-order parameter) for gradient.
        y2 : float, optional
           Y**2 coefficient (2nd-order parameter) for gradient.

        Returns
        -------
        None, or 2-d numpy.array of width/height matching input bbox, containing computed gradient values.
        """

        gradient = None
        if b is not None:  # Don't fit for other gradient parameters if the intercept is not allowed.
            y, x = in_x[0, :], in_x[1, :]
            gradient = np.full_like(x, b, dtype='float64')
            if x1 is not None:
                gradient += x1 * x
            if y1 is not None:
                gradient += y1 * y
            if xy is not None:
                gradient += xy * (x * y)
            if x2 is not None:
                gradient += x2 * (x * x)
            if y2 is not None:
                gradient += y2 * (y * y)
        return gradient

    @staticmethod
    def genStarModel(bbox, psf, xcen, ycen, flux):
        """Generate model (2-d array) of a 'star' (single PSF) centered at given coordinates

        Parameters
        ----------
        bbox : afw.geom.BoundingBox
           Bounding box marking pixel coordinates for generated model
        psf : afw.detection.Psf
           Psf model used to generate the 'star'
        xcen : float
           Desired x-centroid of the 'star'
        ycen : float
           Desired y-centroid of the 'star'
        flux : float
           Desired flux of the 'star'

        Returns
        -------
        2-d numpy.array of width/height matching input bbox, containing PSF with given centroid and flux
        """

        # Generate the psf image, normalize to flux
        psf_img = psf.computeImage(Point2D(xcen, ycen)).convertF()
        psf_img_sum = np.nansum(psf_img.getArray())
        psf_img *= (flux/psf_img_sum)

        # Clip the PSF image bounding box to fall within the footprint bounding box
        psf_box = psf_img.getBBox()
        psf_box.clip(bbox)
        psf_img = ImageF(psf_img, psf_box, PARENT)

        # Then actually crop the psf image.
        # Usually not necessary, but if the dipole is near the edge of the image...
        # Would be nice if we could compare original pos_box with clipped pos_box and
        #     see if it actually was clipped.
        p_Im = ImageF(bbox)
        tmpSubim = ImageF(p_Im, psf_box, PARENT)
        tmpSubim += psf_img

        return p_Im

    @staticmethod
    def genDipoleModel(x, flux, xcenPos, ycenPos, xcenNeg, ycenNeg, fluxNeg=None,
                       b=None, x1=None, y1=None, xy=None, x2=None, y2=None,
                       bNeg=None, x1Neg=None, y1Neg=None, xyNeg=None, x2Neg=None, y2Neg=None,
                       debug=False, **kwargs):
        """
        Generate dipole model with given parameters.
        This is the functor whose sum-of-squared difference from data is minimized by `lmfit`.

        Parameters
        ----------
        x : numpy.array
           Input independent variable. Used here as the grid on which to compute the background
           gradient model.
        flux : float
           Desired flux of the positive lobe of the dipole
        xcenPos : float
           Desired x-centroid of the positive lobe of the dipole
        ycenPos : float
           Desired y-centroid of the positive lobe of the dipole
        xcenNeg : float
           Desired x-centroid of the negative lobe of the dipole
        ycenNeg : float
           Desired y-centroid of the negative lobe of the dipole
        fluxNeg : float, optional
           Desired flux of the negative lobe of the dipole, set to 'flux' if None
        b, x1, y1, xy, x2, y2 : float, optional
           Gradient parameters for positive lobe, see `genBgGradientModel`.
        bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg : float, optional
           Gradient parameters for negative lobe, see `genBgGradientModel`. Set to the corresponding
           positive values if None.
        **kwargs :
           Keyword arguments passed through `lmfit` and used by this functor:
           {
              psf : afw.detection.Psf
                 Psf model used to generate the 'star'
              rel_weight : float
                 Used to signify if positive/negative images are to be included (!= 0. if yes)
              bbox : afw.geom.BoundingBox
                 Bounding box containing region to be modelled
           }

        Returns
        -------
        numpy.array of width/height matching input bbox, containing dipole model with
        given centroids and flux(es). If rel_weight = 0., this is a 2-d array with dimensions
        matching those of bbox; otherwise a stack of three such arrays, representing the dipole
        (diffim), positive and negative images respectively.
        """

        psf = kwargs.get('psf')
        rel_weight = kwargs.get('rel_weight')  # if > 0, we're including pre-sub. images
        fp = kwargs.get('footprint')
        bbox = fp.getBBox()

        if fluxNeg is None:
            fluxNeg = flux

        if debug:
            print('%.2f %.2f %.2f %.2f %.2f %.2f' % (flux, fluxNeg, xcenPos, ycenPos, xcenNeg, ycenNeg))
            if x1 is not None:
                print('     %.2f %.2f %.2f' % (b, x1, y1))
            if xy is not None:
                print('     %.2f %.2f %.2f' % (xy, x2, y2))

        posIm = DipoleFitAlgorithm.genStarModel(bbox, psf, xcenPos, ycenPos, flux)
        negIm = DipoleFitAlgorithm.genStarModel(bbox, psf, xcenNeg, ycenNeg, fluxNeg)

        in_x = x
        if in_x is None:  # use the footprint to generate the input grid
            y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
            in_x = np.array([x, y]) * 1.

        gradient = DipoleFitAlgorithm.genBgGradientModel(in_x, b, x1, y1, xy, x2, y2)
        gradientNeg = gradient
        if bNeg is not None:
            gradientNeg = DipoleFitAlgorithm.genBgGradientModel(in_x, bNeg, x1Neg, y1Neg, xyNeg, x2Neg, y2Neg)

        if gradient is not None:
            posIm.getArray()[:, :] += gradient
            negIm.getArray()[:, :] += gradientNeg

        # Generate the diffIm model
        diffIm = ImageF(bbox)
        diffIm += posIm
        diffIm -= negIm

        zout = diffIm.getArray()
        if rel_weight > 0.:
            zout = np.append([zout], [posIm.getArray(), negIm.getArray()], axis=0)

        return zout

    def fitDipoleImpl(self, source, tol=1e-7, rel_weight=0.5,
                      fitBgGradient=True, bgGradientOrder=1, centroidRangeInSigma=5.,
                      separateNegParams=True, verbose=False, display=False):
        """
        Fit a dipole model to an input difference image (actually,
        subimage bounded by the input source's footprint) and
        optionally constrain the fit using the pre-subtraction images
        posImage and negImage.

        Parameters
        ----------
        See `fitDipole()`

        Returns
        -------
        `lmfit.MinimizerResult` object containing the fit parameters and other information.

        """

        fp = source.getFootprint()
        bbox = fp.getBBox()
        subim = MaskedImageF(self.diffim.getMaskedImage(), bbox, PARENT)

        z = diArr = subim.getArrays()[0]
        weights = 1. / subim.getArrays()[2]  # get the weights (=1/variance)
        if self.posImage is not None and rel_weight > 0.:
            posSubim = MaskedImageF(self.posImage.getMaskedImage(), bbox, PARENT)
            negSubim = MaskedImageF(self.negImage.getMaskedImage(), bbox, PARENT)
            z = np.append([z], [posSubim.getArrays()[0],
                                negSubim.getArrays()[0]], axis=0)
            # Weight the pos/neg images by rel_weight relative to the diffim
            weights = np.append([weights], [1. / posSubim.getArrays()[2] * rel_weight,
                                            1. / negSubim.getArrays()[2] * rel_weight], axis=0)

        psfSigma = self.diffim.getPsf().computeShape().getDeterminantRadius()

        # Create the lmfit model (lmfit uses scipy 'leastsq' option by default - Levenberg-Marquardt)
        gmod = DipoleFitAlgorithm.lmfit.Model(DipoleFitAlgorithm.genDipoleModel, verbose=verbose)

        # Add the constraints for centroids, fluxes.
        # starting constraint - near centroid of footprint
        fpCentroid = np.array([fp.getCentroid().getX(), fp.getCentroid().getY()])
        cenNeg = cenPos = fpCentroid

        pks = fp.getPeaks()
        if len(pks) >= 1:
            cenPos = pks[0].getF()    # if individual (merged) peaks were detected, use those
        if len(pks) >= 2:
            cenNeg = pks[1].getF()

        # For close/faint dipoles the starting locs (min/max) might be way off, let's help them a bit.
        # First assume dipole is not separated by more than 5*psfSigma.
        centroidRange = psfSigma * centroidRangeInSigma

        # Note - this may be a cheat to assume the dipole is centered in center of the footprint.
        if np.sum(np.sqrt((np.array(cenPos) - fpCentroid)**2.)) > centroidRange:
            cenPos = fpCentroid
        if np.sum(np.sqrt((np.array(cenNeg) - fpCentroid)**2.)) > centroidRange:
            cenPos = fpCentroid

        # parameter hints/constraints: https://lmfit.github.io/lmfit-py/model.html#model-param-hints-section
        # might make sense to not use bounds -- see http://lmfit.github.io/lmfit-py/bounds.html
        # also see this discussion -- https://github.com/scipy/scipy/issues/3129
        gmod.set_param_hint('xcenPos', value=cenPos[0],
                            min=cenPos[0]-centroidRange, max=cenPos[0]+centroidRange)
        gmod.set_param_hint('ycenPos', value=cenPos[1],
                            min=cenPos[1]-centroidRange, max=cenPos[1]+centroidRange)
        gmod.set_param_hint('xcenNeg', value=cenNeg[0],
                            min=cenNeg[0]-centroidRange, max=cenNeg[0]+centroidRange)
        gmod.set_param_hint('ycenNeg', value=cenNeg[1],
                            min=cenNeg[1]-centroidRange, max=cenNeg[1]+centroidRange)

        # Estimate starting flux. This strongly affects runtime performance so we want to make it close.
        # I tried many possibilities, the best one seems to be just using sum(abs(diffIm))/2.
        # I am leaving the other tries here (commented out) for posterity and possible future improvements...

        # Value to convert peak value to total flux based on flux within psf
        # psfImg = diffim.getPsf().computeImage()
        # pkToFlux = np.nansum(psfImg.getArray()) / diffim.getPsf().computePeak()

        # bg = np.nanmedian(diArr) # Compute the dipole background (probably very close to zero)
        # startingPk = np.nanmax(diArr) - bg   # use the dipole peak for an estimate.
        # posFlux, negFlux = startingPk * pkToFlux, -startingPk * pkToFlux

        # if len(pks) >= 1:
        #     posFlux = pks[0].getPeakValue() * pkToFlux
        # if len(pks) >= 2:
        #     negFlux = pks[1].getPeakValue() * pkToFlux

        # Use the (flux under the dipole)*5 for an estimate.
        # Lots of testing showed that having startingFlux be too high was better than too low.
        startingFlux = np.nansum(np.abs(diArr) - np.nanmedian(np.abs(diArr))) * 5.
        posFlux = negFlux = startingFlux

        # TBD: set max. flux limit?
        gmod.set_param_hint('flux', value=posFlux, min=0.1)

        if separateNegParams:
            # TBD: set max negative lobe flux limit?
            gmod.set_param_hint('fluxNeg', value=np.abs(negFlux), min=0.1)

        # Fixed parameters (dont fit for them if there are no pre-sub images or no gradient fit requested):
        if (rel_weight > 0. and fitBgGradient):
            if bgGradientOrder >= 0:
                gmod.set_param_hint('b', value=0.)
                if separateNegParams:
                    gmod.set_param_hint('bNeg', value=0.)
            if bgGradientOrder >= 1:
                gmod.set_param_hint('x1', value=0.)
                gmod.set_param_hint('y1', value=0.)
                if separateNegParams:
                    gmod.set_param_hint('x1Neg', value=0.)
                    gmod.set_param_hint('y1Neg', value=0.)
            if bgGradientOrder >= 2:
                gmod.set_param_hint('xy', value=0.)
                gmod.set_param_hint('x2', value=0.)
                gmod.set_param_hint('y2', value=0.)
                if separateNegParams:
                    gmod.set_param_hint('xyNeg', value=0.)
                    gmod.set_param_hint('x2Neg', value=0.)
                    gmod.set_param_hint('y2Neg', value=0.)

        y, x = np.mgrid[bbox.getBeginY():bbox.getEndY(), bbox.getBeginX():bbox.getEndX()]
        in_x = np.array([x, y]).astype(np.float)

        # I'm not sure about the variance planes in the diffim (or convolved pre-sub. images
        # for that matter) so for now, let's just do an un-weighted least-squares fit
        # (override weights computed above).
        weights = 1.
        if self.posImage is not None and rel_weight > 0.:
            weights = np.array([np.ones_like(diArr), np.ones_like(diArr)*rel_weight,
                                np.ones_like(diArr)*rel_weight])

        # Note that although we can, we're not required to set initial values for params here,
        # since we set their param_hint's above.
        # add "method" param to not use 'leastsq' (==levenberg-marquardt), e.g. "method='nelder'"
        result = gmod.fit(z, weights=weights, x=in_x,
                          verbose=verbose,
                          fit_kws={'ftol': tol, 'xtol': tol, 'gtol': tol, 'maxfev': 250},  # see scipy docs
                          psf=self.diffim.getPsf(),  # hereon: kwargs that get passed to genDipoleModel()
                          rel_weight=rel_weight,
                          footprint=fp)

        # Probably never wanted - also this takes a long time (longer than the fit!)
        # This is how to get confidence intervals out:
        #    https://lmfit.github.io/lmfit-py/confidence.html and
        #    http://cars9.uchicago.edu/software/python/lmfit/model.html
        if verbose:  # the ci_report() seems to fail if neg params are constrained -- TBD why.
            print(result.fit_report(show_correl=False))
            if separateNegParams:
                print(result.ci_report())

        return result

    def fitDipole(self, source, tol=1e-7, rel_weight=0.1,
                  fitBgGradient=True, centroidRangeInSigma=5., separateNegParams=True,
                  bgGradientOrder=1, verbose=False, display=False, return_fitObj=False):
        """
        Wrapper around `fitDipoleImpl()` which performs the fit of a dipole
        model to an input difference image (actually, subimage bounded
        by the input source's footprint) and optionally constrain the
        fit using the pre-subtraction images posImage and
        negImage. Wraps the output into a `resultsOutput` object after
        computing additional statistics such as orientation and SNR.

        Parameters
        ----------
        source : afw.table.SourceRecord
           Record containing the (merged) dipole source footprint detected on the diffim
        tol : float, optional
           Tolerance parameter for scipy.leastsq() optimization
        rel_weight : float, optional
           Weighting of posImage/negImage relative to the diffim in the fit
        fitBgGradient : bool, optional
           Fit linear background gradient in posImage/negImage?
        bgGradientOrder : int, optional
           Desired polynomial order of background gradient (allowed are [0,1,2])
        centroidRangeInSigma : float, optional
           Allowed window of centroid parameters relative to peak in input source footprint
        separateNegParams : bool, optional
           Fit separate parameters to the flux and background gradient in the negative images?
        verbose : bool, optional
           Be verbose
        display : bool, optional
           Display input data, best fit model(s) and residuals in a matplotlib window.
        return_fitObj : bool, optional
           In addition to the `resultsOutput` object, also returns the
           `lmfit.MinimizerResult` object for debugging.

        Returns
        -------
        `resultsOutput` object containing the fit parameters and other information.
        `lmfit.MinimizerResult` object if `return_fitObj` is True.
        """

        fitResult = self.fitDipoleImpl(
            source, tol=tol, rel_weight=rel_weight, fitBgGradient=fitBgGradient,
            centroidRangeInSigma=centroidRangeInSigma, separateNegParams=separateNegParams,
            bgGradientOrder=bgGradientOrder, verbose=verbose, display=display)

        fitParams = fitResult.best_values
        if fitParams['flux'] <= 1.:   # usually around 0.1 -- the minimun flux allowed -- i.e. bad fit.
            out = DipoleFitAlgorithm.resultsOutput(
                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan)
            if return_fitObj:  # for debugging
                return out, fitResult
            return out

        centroid = ((fitParams['xcenPos'] + fitParams['xcenNeg']) / 2.,
                    (fitParams['ycenPos'] + fitParams['ycenNeg']) / 2.)
        dx, dy = fitParams['xcenPos'] - fitParams['xcenNeg'], fitParams['ycenPos'] - fitParams['ycenNeg']
        angle = np.arctan2(dy, dx) / np.pi * 180.   # convert to degrees (should keep as rad?)

        # Exctract flux value, compute signalToNoise from flux/variance_within_footprint
        # Also extract the stderr of flux estimate.
        def computeSumVariance(exposure, footprint):
            box = footprint.getBBox()
            subim = MaskedImageF(exposure.getMaskedImage(), box, PARENT)
            return np.sqrt(np.nansum(subim.getArrays()[1][:, :]))

        fluxVal = fluxVar = fitParams['flux']
        fluxErr = fluxErrNeg = fitResult.params['flux'].stderr
        if self.posImage is not None:
            fluxVar = computeSumVariance(self.posImage, source.getFootprint())
        else:
            fluxVar = computeSumVariance(self.diffim, source.getFootprint())

        fluxValNeg, fluxVarNeg = fluxVal, fluxVar
        if separateNegParams:
            fluxValNeg = fitParams['fluxNeg']
            fluxErrNeg = fitResult.params['fluxNeg'].stderr
        if self.negImage is not None:
            fluxVarNeg = computeSumVariance(self.negImage, source.getFootprint())

        try:
            signalToNoise = np.sqrt((fluxVal/fluxVar)**2 + (fluxValNeg/fluxVarNeg)**2)
        except:  # catch divide by zero - should never happen.
            signalToNoise = np.nan

        out = DipoleFitAlgorithm.resultsOutput(
            fitParams['xcenPos'], fitParams['ycenPos'], fitParams['xcenNeg'], fitParams['ycenNeg'],
            fluxVal, -fluxValNeg, fluxErr, fluxErrNeg, centroid[0], centroid[1], angle,
            signalToNoise, fitResult.chisqr, fitResult.redchi)

        if return_fitObj:  # for debugging
            return out, fitResult
        return out
