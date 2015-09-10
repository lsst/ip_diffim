#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy

import lsst.pex.config as pexConfig
from .dipoleMeasurement import DipoleAnalysis
from lsst.meas.base.pluginRegistry import register
from lsst.meas.base.sfm import SingleFramePluginConfig, SingleFramePlugin

__all__ = ["SingleFrameClassificationDipoleConfig", "SingleFrameClassificationDipolePlugin"]


class SingleFrameClassificationDipoleConfig(SingleFramePluginConfig):
    """Configuration for classification of detected diaSources as dipole or not"""
    minSn = pexConfig.Field(
        doc="Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole",
        dtype=float, default=numpy.sqrt(2) * 5.0,
        )
    maxFluxRatio = pexConfig.Field(
        doc = "Maximum flux ratio in either lobe to be considered a dipole",
        dtype = float, default = 0.65
        )


@register("ip_diffim_ClassificationDipole")
class SingleFrameClassificationDipolePlugin(SingleFramePlugin):
    """
    A binary measure of the whether a diasource is a dipole.
    """

    ConfigClass = SingleFrameClassificationDipoleConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CLASSIFY_ORDER

    def __init__(self, config, name, schema, metadata):
        SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.dipoleAnalysis = DipoleAnalysis()
        self._ClassificationFlag = "classification_dipole"
        self.keyProbability = schema.addField(name + "_value", type="D",
                                              doc="Set to 1 for dipoles, else 0.")
        self.keyFlag = schema.addField(name + "_flag", type="Flag", doc="Set to 1 for any fatal failure.")

    def measure(self, measRecord, exposure):
            passesSn = self.dipoleAnalysis.getSn(measRecord) > self.config.minSn
            negFlux = numpy.abs(measRecord.get("ip_diffim_PsfDipoleFlux_neg_flux"))
            negFluxFlag =  measRecord.get("ip_diffim_PsfDipoleFlux_neg_flag")
            posFlux = numpy.abs(measRecord.get("ip_diffim_PsfDipoleFlux_pos_flux"))
            posFluxFlag = measRecord.get("ip_diffim_PsfDipoleFlux_pos_flag")
            if numpy.isnan(negFlux) or numpy.isnan(posFlux) or negFluxFlag or posFluxFlag:
                self.fail(measRecord)
            totalFlux = negFlux + posFlux
            passesFluxNeg = (negFlux / totalFlux) < self.config.maxFluxRatio
            passesFluxPos = (posFlux / totalFlux) < self.config.maxFluxRatio
            if (passesSn and passesFluxPos and passesFluxNeg):
                val = 1.0
            else:
                val = 0.0
            measRecord.set(self.keyProbability, val)

    def fail(self, measRecord, error=None):
        # Override fail() to do nothing in the case of an exception.  We should be setting a flag
        # instead.
        measRecord.set(self.keyFlag, True)