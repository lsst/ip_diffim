.. lsst-task-topic:: lsst.ip.diffim.DetectAndMeasureTask

#####################
DetectAndMeasureTask
#####################

DetectAndMeasureTask performs source detection and execution of measurement plugins in difference images.

.. _lsst.ip.diffim.DetectAndMeasureTask-description:

Description
===========

This task cross-convolves the difference image with a PSF approximation to obtain maximum-likelihood source candidates.
Measurements are run on all the detections, including positive and negative flux candidates.
Close sources will be combined, enabling trailed sources as well as dipole-like source detection.

The plugins used by this task by default are meant to provide basic measurements for all types of transient sources: point sources, trailed sources and dipoles as well.
Additional forced measurements are provided through PSF flux photometry measurements and centroid determination plugins.
Finally, the plugins make use of the image mask planes to provide a wide range of source flags, including streak flagging and injected source flagging.

.. _lsst.ip.diffim.DetectAndMeasureTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.DetectAndMeasureTask

.. _lsst.ip.diffim.DetectAndMeasureTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.DetectAndMeasureTask

.. _lsst.ip.diffim.DetectAndMeasureTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.DetectAndMeasureTask

.. _lsst.ip.diffim.DetectAndMeasureTask-debug:

.. Debugging
.. =========


.. The ``pipetask`` command line interface supports a ``--debug`` flag to import
.. ``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
.. files.
.. The available variables in DetectAndMeasureTask include:


.. display : `bool`
..     Enable debug display output.
.. maskTransparency : `float`
..     Transparency of mask planes in the output display.
.. displayDiaSources : `bool`
..     Show exposure with dipole fitting results.
