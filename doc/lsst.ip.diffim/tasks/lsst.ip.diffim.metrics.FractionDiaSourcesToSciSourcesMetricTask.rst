.. lsst-task-topic:: lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask

########################################
FractionDiaSourcesToSciSourcesMetricTask
########################################

``FractionDiaSourcesToSciSourcesMetricTask`` computes the fraction of science sources that become DIASources when processed through image differencing (as the ``ip_diffim.fracDiaSourcesToSciSources`` metric).
It requires source catalogs as input, and can operate at image-level or coarser granularity.

.. _lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask-summary:

Processing summary
==================

``FractionDiaSourcesToSciSourcesMetricTask`` reads source catalogs (``src`` and ``deepDiff_diaSrc`` datasets, by default) and counts the number of sources.
It then computes the ratio of DIASources to science sources for those catalogs.

.. _lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask

.. _lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask-butler:

Butler datasets
===============

Input datasets
--------------

:lsst-config-field:`~_lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricConfig.sciSources`
    The catalog type for science sources (default: ``src``).

:lsst-config-field:`~_lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricConfig.diaSources`
    The catalog type for DIASources (default: ``deepDiff_diaSrc``).

.. _lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask

.. _lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.metrics.FractionDiaSourcesToSciSourcesMetricTask
