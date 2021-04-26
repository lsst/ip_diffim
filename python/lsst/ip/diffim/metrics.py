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


import numpy as np
import astropy.units as u

from lsst.pipe.base import Struct, connectionTypes
from lsst.verify import Measurement
from lsst.verify.gen2tasks import register
from lsst.verify.tasks import MetricTask, MetricConfig, MetricConnections, \
    MetricComputationError


class NumberSciSourcesMetricConnections(
        MetricConnections,
        defaultTemplates={"package": "ip_diffim",
                          "metric": "numSciSources"},
        dimensions={"instrument", "visit", "detector"},
):
    sources = connectionTypes.Input(
        doc="The catalog of science sources.",
        name="src",
        storageClass="SourceCatalog",
        dimensions={"instrument", "visit", "detector"},
    )


class NumberSciSourcesMetricConfig(
        MetricConfig,
        pipelineConnections=NumberSciSourcesMetricConnections):
    pass


@register("numSciSources")
class NumberSciSourcesMetricTask(MetricTask):
    """Task that computes the number of cataloged non-sky science sources.

    Notes
    -----
    The task excludes any sky sources in the catalog, but it does not require
    that the catalog include a ``sky_sources`` column.
    """
    _DefaultName = "numSciSources"
    ConfigClass = NumberSciSourcesMetricConfig

    def run(self, sources):
        """Count the number of non-sky science sources.

        Parameters
        ----------
        sources : `lsst.afw.table.SourceCatalog` or `None`
            A science source catalog, which may be empty or `None`.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A `~lsst.pipe.base.Struct` containing the following component:

            ``measurement``
                the total number of non-sky science sources (`lsst.verify.Measurement`
                or `None`)
        """
        if sources is not None:
            nSciSources = _countRealSources(sources)
            meas = Measurement(self.config.metricName, nSciSources * u.count)
        else:
            self.log.info("Nothing to do: no catalogs found.")
            meas = None
        return Struct(measurement=meas)


class FractionDiaSourcesToSciSourcesMetricConnections(
        MetricTask.ConfigClass.ConnectionsClass,
        dimensions={"instrument", "visit", "detector"},
        defaultTemplates={"coaddName": "deep",
                          "fakesType": "",
                          "package": "ip_diffim",
                          "metric": "fracDiaSourcesToSciSources"}):
    sciSources = connectionTypes.Input(
        doc="The catalog of science sources.",
        name="src",
        storageClass="SourceCatalog",
        dimensions={"instrument", "visit", "detector"},
    )
    diaSources = connectionTypes.Input(
        doc="The catalog of DIASources.",
        name="{fakesType}{coaddName}Diff_diaSrc",
        storageClass="SourceCatalog",
        dimensions={"instrument", "visit", "detector"},
    )


class FractionDiaSourcesToSciSourcesMetricConfig(
        MetricTask.ConfigClass,
        pipelineConnections=FractionDiaSourcesToSciSourcesMetricConnections):
    pass


@register("fracDiaSourcesToSciSources")
class FractionDiaSourcesToSciSourcesMetricTask(MetricTask):
    """Task that computes the ratio of difference image sources to science
    sources in an image, visit, etc.

    Notes
    -----
    The task excludes any sky sources in the direct source catalog, but it
    does not require that either catalog include a ``sky_sources`` column.
    """
    _DefaultName = "fracDiaSourcesToSciSources"
    ConfigClass = FractionDiaSourcesToSciSourcesMetricConfig

    def run(self, sciSources, diaSources):
        """Compute the ratio of DIASources to non-sky science sources.

        Parameters
        ----------
        sciSources : `lsst.afw.table.SourceCatalog` or `None`
            A science source catalog, which may be empty or `None`.
        diaSources : `lsst.afw.table.SourceCatalog` or `None`
            A DIASource catalog for the same unit of processing
            as ``sciSources``.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            A `~lsst.pipe.base.Struct` containing the following component:

            ``measurement``
                the ratio (`lsst.verify.Measurement` or `None`)
        """
        if diaSources is not None and sciSources is not None:
            nSciSources = _countRealSources(sciSources)
            nDiaSources = _countRealSources(diaSources)
            metricName = self.config.metricName
            if nSciSources <= 0:
                raise MetricComputationError(
                    "No science sources found; ratio of DIASources to science sources ill-defined.")
            else:
                meas = Measurement(metricName, nDiaSources / nSciSources * u.dimensionless_unscaled)
        else:
            self.log.info("Nothing to do: no catalogs found.")
            meas = None
        return Struct(measurement=meas)


def _countRealSources(catalog):
    """Return the number of valid sources in a catalog.

    At present, this definition excludes sky sources. If a catalog does not
    have a ``sky_source`` flag, all sources are assumed to be non-sky.

    Parameters
    ----------
    `catalog` : `lsst.afw.table.SourceCatalog`
        The catalog of sources to count.

    Returns
    -------
    count : `int`
        The number of sources that satisfy the criteria.
    """
    if "sky_source" in catalog.schema:
        return np.count_nonzero(catalog["sky_source"] == False)  # noqa: E712
    else:
        return len(catalog)
