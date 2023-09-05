.. py:currentmodule:: lsst.ip.diffim

.. _lsst.ip.diffim:

##############
lsst.ip.diffim
##############

.. Paragraph that describes what this Python module does and links to related modules and frameworks.

The ``lsst.ip.diffim`` module provides algorithms and utilities for astronomical **image differencing** and **transient detection**.

.. _lsst.ip.diffimg-overview:

Overview of lsst.ip.diffim
==========================

.. toctree linking to the basic overview of the package

.. toctree::
   :maxdepth: 1 

.. toctree::
   :maxdepth: 1 

   overview_ipdiffim
   AL_implementation


.. _lsst.ip.diffim-using:

Using lsst.ip.diffim
====================

.. toctree linking to topics related to using the module's APIs.

.. toctree
   :maxdepth: 1

.. toctree
   :maxdepth: 1

Production pipelines
--------------------

This package is used in production pipelines like ``lsst.ap.pipe`` and ``lsst.drp.pipe``, designed for processing large collections of survey data.

For an example of how the two primary `~lsst.pipe.base.PipelineTasks`  in this package (`~lsst.ip.diffim.AlardLuptonSubtractTask` and `~lsst.ip.diffim.DetectAndMeasureTask`) are used in a production pipeline, see `ap_pipe/pipelines/ApPipe.yaml`.


Single frame processing
-----------------------

Single image frame usage can be carried out by importing the `~lsst.ip.diffim` as a library of tasks, as the following example shows:

.. code-block:: python

   from lsst.ip.diffim import subtractImages

   ...
   subtract_task = subtractImages.AlardLuptonSubtractTask()

   subtraction = subtract_task.run(
      warped_template, science, science_source_catalog)

For source detection a similar pattern works:

.. code-block:: python

   from lsst.ip.diffim import DetectAndMeasureTask

   ...
   detect_task = DetectAndMeasureTask()

   detect_and_measure = detect_and_measure_task.run(
      science, subtraction.matchedTemplate, subtraction.difference)


.. _lsst.ip.diffim-contributing:

Contributing
============

``lsst.ip.diffim`` is developed at https://github.com/lsst/ip_diffim.
You can find Jira issues for this module under the `ip_diffim <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20component%20%3D%20ip_diffim>`_ component.

.. _lsst.ip.diffim-command-line-taskref:

Task reference
==============

.. _lsst.ip.diffim-pipeline-tasks:

Pipeline tasks
--------------

.. lsst-pipelinetasks::
   :root: lsst.ip.diffim

.. _lsst.ip.diffim-tasks:

Tasks
-----

.. lsst-tasks::
   :root: lsst.ip.diffim
   :toctree: tasks

.. _lsst.ip.diffim-configs:

Configurations
--------------

.. lsst-configs::
   :root: lsst.ip.diffim
   :toctree: config

.. _lsst.ip.diffim-pyapi:

Python API reference
====================

.. automodapi:: lsst.ip.diffim
   :no-main-docstr:
   :skip: wrapSimpleAlgorithm

.. automodapi:: lsst.ip.diffim.metrics
   :no-main-docstr:
   :no-inheritance-diagram:
