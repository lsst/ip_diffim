.. lsst-task-topic:: lsst.ip.diffim.DipoleMeasurementTask

#####################
DipoleMeasurementTask
#####################

DipoleMeasurementTask measures sources on images to identify likely dipoles.
This task is purely a configuration wrapper around `~lsst.meas.base.SingleFrameMeasurementTask`, with helpers to performe dipole-specific measurements.

.. _lsst.ip.diffim.DipoleMeasurementTask-description:

Description
===========

The plugins used by this task by default allow the user to test the hypothesis that the
Source is a dipole. This includes a set of measurements derived from
intermediate base classes DipoleCentroidAlgorithm and DipoleFluxAlgorithm.
Their respective algorithm control classes are defined in
DipoleCentroidControl and DipoleFluxControl. Each centroid and flux
measurement will have _neg (negative) and _pos (positive lobe) fields.

The first set of measurements uses the SDSS algorithm for centroid and flux
measurements.
This centroid will be used as a fallback, if the below dipole fit fails or if the source does not have a positive and negative peak to attempt a dipole fit to.

The second set of measurements undertakes a joint-Psf model on the negative and
positive lobe simultaneously. This fit simultaneously solves for the negative
and positive lobe centroids and fluxes using non-linear least squares
minimization. The fields are stored in table elements
ip_diffim_PsfDipoleFlux*.

.. _lsst.ip.diffim.DipoleMeasurementTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.DipoleMeasurementTask

.. _lsst.ip.diffim.DipoleMeasurementTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.DipoleMeasurementTask

.. _lsst.ip.diffim.DipoleMeasurementTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.DipoleMeasurementTask

.. _lsst.ip.diffim.DipoleMeasurementTask-debug:

Debugging
=========


The ``pipetask`` command line interface supports a ``--debug`` flag to import
``debug.py`` from your PYTHONPATH; see :ref:`lsstDebug` for more about ``debug.py``
files.
The available variables in DipoleMeasurementTask include:


display : `bool`
    Enable debug display output.
maskTransparency : `float`
    Transparency of mask planes in the output display.
displayDiaSources : `bool`
    Show exposure with dipole fitting results.
