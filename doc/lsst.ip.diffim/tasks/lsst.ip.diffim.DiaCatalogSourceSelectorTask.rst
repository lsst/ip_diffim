.. lsst-task-topic:: lsst.ip.diffim.DiaCatalogSourceSelectorTask

############################
DiaCatalogSourceSelectorTask
############################

A task that selects sources for Kernel candidates.

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-description:

Description
===========

A naive star selector based on second moments. Use with caution.

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.DiaCatalogSourceSelectorTask

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.DiaCatalogSourceSelectorTask

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.DiaCatalogSourceSelectorTask

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-debug:

Debugging
=========

The ``pipetask`` command line interface supports a ``--debug`` flag to import
``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
files.
The available variables in DiaCatalogSourceSelectorTask include:

display : `bool`
    if True display debug information
displayExposure : `bool`
    if True display exposure
pauseAtEnd `bool`
    if True wait after displaying everything and wait for user input
