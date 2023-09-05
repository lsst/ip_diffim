.. lsst-task-topic:: lsst.ip.diffim.AlardLuptonSubtractTask

#######################
AlardLuptonSubtractTask
#######################

AlardLuptonSubtractTask performs image differencing, including decorrelation.


.. _lsst.ip.diffim.AlardLuptonSubtractTask-description:

Description
===========

This task takes the pair of science and template exposures, together with the science image source catalog, and produces difference and matched template exposures as outputs.

The process involves selecting sources from the science exposure source catalog, obtaining the convolution kernel, applying this to the `target` image, and subtracting the images.

As discussed in the general package overview the task has two possible modes: one convolving the `template` image and the alternative convolving the `science` image.
Both modes share the preparation sub-tasks, but the method that performs the convolution is not the same.
The parameter that controls this behavior is ``mode`` (check `~lsst.ip.diffim.AlardLuptonSubtractConfig`), and it allows ``convolveScience``, ``convolveTemplate`` and finally ``auto`` which means that the decision is made at runtime automatically.

Finally, mask planes are set on the difference and matched template outputs based on the masks of the input science and template images.


.. _lsst.ip.diffim.AlardLuptonSubtractTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.AlardLuptonSubtractTask

.. _lsst.ip.diffim.AlardLuptonSubtractTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.AlardLuptonSubtractTask

.. _lsst.ip.diffim.AlardLuptonSubtractTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.AlardLuptonSubtractTask

.. _lsst.ip.diffim.AlardLuptonSubtractTask-debug:

.. Debugging
.. =========


.. The ``pipetask`` command line interface supports a ``--debug`` flag to import
.. ``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
.. files.
.. The available variables in AlardLuptonSubtractTask include:


.. display : `bool`
..     Enable debug display output.
.. maskTransparency : `float`
..     Transparency of mask planes in the output display.
.. displayDiaSources : `bool`
..     Show exposure with dipole fitting results.
