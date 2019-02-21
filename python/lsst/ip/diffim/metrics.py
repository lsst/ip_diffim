# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__all__ = [
    "NumberSciSourcesMetricTask", "NumberSciSourcesMetricConfig",
    "FractionDiaSourcesToSciSourcesMetricTask", "FractionDiaSourcesToSciSourcesMetricConfig",
]


import astropy.units as u

from lsst.pipe.base import Struct, InputDatasetField
from lsst.verify import Measurement
from lsst.verify.gen2tasks import MetricTask, register
from lsst.verify.tasks import MetricComputationError


class NumberSciSourcesMetricConfig(MetricTask.ConfigClass):
    sources = InputDatasetField(
        doc="The catalog of science sources.",
        name="src",
        storageClass="SourceCatalog",
        dimensions={"Instrument", "Exposure", "Detector"},
    )


@register("numSciSources")
class NumberSciSourcesMetricTask(MetricTask):
    """Task that computes the number of cataloged science sources.
    """
    _DefaultName = "numSciSources"
    ConfigClass = NumberSciSourcesMetricConfig

    def run(self, sources):
        """Count the number of science sources.

        Parameters
        ----------
        sources : iterable of `lsst.afw.table.SourceCatalog`
            A collection of science source catalogs, one for each unit of
            processing to be incorporated into this metric. Its elements may
            be `None` to represent missing data.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A `~lsst.pipe.base.Struct` containing the following component:

            ``measurement``
                the total number of science sources (`lsst.verify.Measurement`
                or `None`)
        """
        nSciSources = 0
        inputData = False
        for catalog in sources:
            if catalog is not None:
                nSciSources += len(catalog)
                inputData = True

        if inputData:
            meas = Measurement(self.getOutputMetricName(self.config), nSciSources * u.count)
        else:
            self.log.info("Nothing to do: no catalogs found.")
            meas = None
        return Struct(measurement=meas)

    @classmethod
    def getOutputMetricName(cls, config):
        return "ip_diffim.numSciSources"


class FractionDiaSourcesToSciSourcesMetricConfig(MetricTask.ConfigClass):
    sciSources = InputDatasetField(
        doc="The catalog of science sources.",
        name="src",
        storageClass="SourceCatalog",
        dimensions={"Instrument", "Exposure", "Detector"},
    )
    diaSources = InputDatasetField(
        doc="The catalog of DIASources.",
        name="deepDiff_diaSrc",
        nameTemplate="{coaddName}Diff_diaSrc",
        storageClass="SourceCatalog",
        dimensions={"Instrument", "Exposure", "Detector"},
    )


@register("fracDiaSourcesToSciSources")
class FractionDiaSourcesToSciSourcesMetricTask(MetricTask):
    """Task that computes the ratio of difference image sources to science
    sources in an image, visit, etc.
    """
    _DefaultName = "fracDiaSourcesToSciSources"
    ConfigClass = FractionDiaSourcesToSciSourcesMetricConfig

    def run(self, sciSources, diaSources):
        """Compute the ratio of DIASources to science sources.

        Parameters
        ----------
        sciSources : iterable of `lsst.afw.table.SourceCatalog`
            A collection of science source catalogs, one for each unit of
            processing to be incorporated into this metric. Its elements may
            be `None` to represent missing data.
        diaSources : iterable of `lsst.afw.table.SourceCatalog`
            A collection of difference imaging catalogs similar to
            ``sciSources``.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A `~lsst.pipe.base.Struct` containing the following component:

            ``measurement``
                the ratio (`lsst.verify.Measurement` or `None`)
        """
        nSciSources = 0
        nDiaSources = 0
        inputData = False

        for sciCatalog, diaCatalog in zip(sciSources, diaSources):
            if diaCatalog is not None and sciCatalog is not None:
                nSciSources += len(sciCatalog)
                nDiaSources += len(diaCatalog)
                inputData = True

        if inputData:
            metricName = self.getOutputMetricName(self.config)
            if nSciSources <= 0.0:
                raise MetricComputationError(
                    "No science sources found; ratio of DIASources to science sources ill-defined.")
                meas = Measurement(metricName, 0.0 * u.dimensionless_unscaled)
            else:
                meas = Measurement(metricName, nDiaSources / nSciSources * u.dimensionless_unscaled)
        else:
            self.log.info("Nothing to do: no catalogs found.")
            meas = None
        return Struct(measurement=meas)

    @classmethod
    def getOutputMetricName(cls, config):
        return "ip_diffim.fracDiaSourcesToSciSources"
