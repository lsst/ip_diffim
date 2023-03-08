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

import numpy as np

import lsst.afw.image
import lsst.afw.math
import lsst.geom
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
from lsst.utils.timer import timeMethod

from lsst.ip.diffim.transiNetInterface import TransiNetInterface

__all__ = ["TransiNetSubtractConfig", "TransiNetSubtractTask"]

_dimensions = ("instrument", "visit", "detector")
_defaultTemplates = {"coaddName": "deep", "fakesType": ""}


class SubtractInputConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=_dimensions,
                               defaultTemplates=_defaultTemplates):
    template = connectionTypes.Input(
        doc="Input warped template to subtract.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_templateExp"
    )
    science = connectionTypes.Input(
        doc="Input science exposure to subtract from.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}calexp"
    )


class SubtractImageOutputConnections(lsst.pipe.base.PipelineTaskConnections,
                                     dimensions=_dimensions,
                                     defaultTemplates=_defaultTemplates):
    difference = connectionTypes.Output(
        doc="Result of subtracting convolved template from science image.",
        dimensions=("instrument", "visit", "detector"),
        storageClass="ExposureF",
        name="{fakesType}{coaddName}Diff_differenceTempExp",
    )


class TransiNetSubtractConnections(SubtractInputConnections, SubtractImageOutputConnections):

    def __init__(self, *, config=None):
        super().__init__(config=config)


class TransiNetSubtractConfig(lsst.pipe.base.PipelineTaskConfig,
                              pipelineConnections=TransiNetSubtractConnections):
    requiredTemplateFraction = lsst.pex.config.Field(
        dtype=float,
        default=0.1,
        doc="Abort task if template covers less than this fraction of pixels."
        " Setting to 0 will always attempt image subtraction."
    )


class TransiNetSubtractTask(lsst.pipe.base.PipelineTask):
    """Compute the image difference of a science and template image using
    the TransiNet algorithm (Sedaghat, Mahabal 2018).
    """
    ConfigClass = TransiNetSubtractConfig
    _DefaultName = "transiNetSubtract"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.transiNetInterface = TransiNetInterface('temp_var', 'neighbor')

    @timeMethod
    def run(self, template, science):
        """ Subtract two images.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.

        Returns
        -------
        results : `lsst.pipe.base.Struct`
            ``difference`` : `lsst.afw.image.ExposureF`
                Result of subtracting template and science.

        Raises
        ------
        RuntimeError
            TBD
        lsst.pipe.base.NoWorkFound
            Raised if fraction of good pixels, defined as not having NO_DATA
            set, is less then the configured requiredTemplateFraction
        """
        self._validateExposures(template, science)
        checkTemplateIsSufficient(template, self.log,
                                  requiredTemplateFraction=self.config.requiredTemplateFraction)

        # put the template on the same photometric scale as the science image
        if False:  # TODO: See if we need this or not -- probably not
            photoCalib = template.getPhotoCalib()
            self.log.info("Applying photometric calibration to template: %f", photoCalib.getCalibrationMean())
            template.maskedImage = photoCalib.calibrateImage(template.maskedImage)

        # Crop the template to the science image's bounding box
        template = template[science.getBBox()]

        difference = self.transiNetInterface.infer(template, science)

        return lsst.pipe.base.Struct(difference=difference)

    @staticmethod
    def _validateExposures(template, science):
        """Check that the WCS of the two Exposures match, and the template bbox
        contains the science bbox.

        Parameters
        ----------
        template : `lsst.afw.image.ExposureF`
            Template exposure, warped to match the science exposure.
        science : `lsst.afw.image.ExposureF`
            Science exposure to subtract from the template.

        Raises
        ------
        AssertionError
            Raised if the WCS of the template is not equal to the science WCS,
            or if the science image is not fully contained in the template
            bounding box.
        """
        assert template.wcs == science.wcs,\
            "Template and science exposure WCS are not identical."
        templateBBox = template.getBBox()
        scienceBBox = science.getBBox()

        assert templateBBox.contains(scienceBBox),\
            "Template bbox does not contain all of the science image."


def checkTemplateIsSufficient(templateExposure, logger, requiredTemplateFraction=0.):
    """Raise NoWorkFound if template coverage < requiredTemplateFraction

    Parameters
    ----------
    templateExposure : `lsst.afw.image.ExposureF`
        The template exposure to check
    logger : `lsst.log.Log`
        Logger for printing output.
    requiredTemplateFraction : `float`, optional
        Fraction of pixels of the science image required to have coverage
        in the template.

    Raises
    ------
    lsst.pipe.base.NoWorkFound
        Raised if fraction of good pixels, defined as not having NO_DATA
        set, is less then the configured requiredTemplateFraction
    """
    # Count the number of pixels with the NO_DATA mask bit set
    # counting NaN pixels is insufficient because pixels without data are often intepolated over)
    pixNoData = np.count_nonzero(templateExposure.mask.array
                                 & templateExposure.mask.getPlaneBitMask('NO_DATA'))
    pixGood = templateExposure.getBBox().getArea() - pixNoData
    logger.info("template has %d good pixels (%.1f%%)", pixGood,
                100*pixGood/templateExposure.getBBox().getArea())

    if pixGood/templateExposure.getBBox().getArea() < requiredTemplateFraction:
        message = ("Insufficient Template Coverage. (%.1f%% < %.1f%%) Not attempting subtraction. "
                   "To force subtraction, set config requiredTemplateFraction=0." % (
                       100*pixGood/templateExposure.getBBox().getArea(),
                       100*requiredTemplateFraction))
        raise lsst.pipe.base.NoWorkFound(message)
