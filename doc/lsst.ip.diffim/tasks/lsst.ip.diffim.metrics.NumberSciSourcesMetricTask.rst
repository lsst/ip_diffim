.. lsst-task-topic:: lsst.ip.diffim.metrics.NumberSciSourcesMetricTask

##########################
NumberSciSourcesMetricTask
##########################

``NumberSciSourcesMetricTask`` computes the number of science sources created when processing data through image differencing (as the ``ip_diffim.numSciSources`` metric).
It requires source catalogs (a ``src`` dataset) as input, and can operate at image-level granularity.

.. _lsst.ip.diffim.metrics.NumberSciSourcesMetricTask-summary:

Processing summary
==================

``NumberSciSourcesToSciSourcesMetricTask`` reads source catalogs (``src`` datasets) and adds up the number of primary sources in those catalogs.

.. _lsst.ip.diffim.metrics.NumberSciSourcesMetricTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.metrics.NumberSciSourcesMetricTask

.. _lsst.ip.diffim.metrics.NumberSciSourcesMetricTask-butler:

Butler datasets
===============

Input datasets
--------------

``sources``
    The catalog type for science sources (default: ``src``).

.. _lsst.ip.diffim.metrics.NumberSciSourcesMetricTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.metrics.NumberSciSourcesMetricTask

.. _lsst.ip.diffim.metrics.NumberSciSourcesMetricTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.metrics.NumberSciSourcesMetricTask
