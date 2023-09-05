.. lsst-task-topic:: lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask

##################################
AlardLuptonPreconvolveSubtractTask
##################################

AlardLuptonPreconvolveSubtractTask is a modification of AlardLuptonSubtractTask behavior, and it performs image differencing between the `template` image and a custom kernel cross-correlated `science` image.

.. _lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask-description:

Description
===========

This task takes the pair of science and template exposures, together with the science image source catalog, and produces difference and matched template exposures as outputs.

The process involves pre-convolving the `science` image with a predefined kernel (tipically its own PSF) before running the nominal set of steps for image subtraction, that is, selecting sources to obtain a matching convolution kernel, applying this to the `template` image, and performing the image differencing.

In this case the subtraction produces a score image directly, and it can be used for maximum-likelihood source detection in a straightforward manner. This is done using `~lsst.ip.diffim.DetectAndMeasureScoreTask`.


.. _lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask

.. _lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask

.. _lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask

.. _lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask-debug:

.. Debugging
.. =========


.. The ``pipetask`` command line interface supports a ``--debug`` flag to import
.. ``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
.. files.
.. The available variables in AlardLuptonPreconvolveSubtractTask include:


.. display : `bool`
..     Enable debug display output.
.. maskTransparency : `float`
..     Transparency of mask planes in the output display.
.. displayDiaSources : `bool`
..     Show exposure with dipole fitting results.
