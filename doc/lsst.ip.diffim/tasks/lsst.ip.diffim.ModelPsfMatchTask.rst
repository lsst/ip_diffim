.. lsst-task-topic:: lsst.ip.diffim.ModelPsfMatchTask

##########################
ModelPsfMatchTask
##########################

ImagePsfMatchTask creates a PSF-matching kernel for two models.

Description
===========

This Task differs from ``ImagePsfMatchTask`` in that it matches two
Psf _models_ by realizing them in an Exposure-sized SpatialCellSet and then
inserting each Psf-image pair into ``KernelCandidates``. Because none of the
pairs of sources that are to be matched should be invalid, all sigma clipping
is turned off in ModelPsfMatchConfig.  And because there is no
tracked _variance_ in the Psf images, the debugging and logging QA info should
be interpreted with caution.

One item of note is that the sizes of Psf models are fixed (e.g. its defined as
a 21x21 matrix).  When the Psf-matching kernel is being solved for, the
Psf "image" is convolved with each kernel basis function, leading to a loss of
information around the borders. This pixel loss will be problematic for the
numerical stability of the kernel solution if the size of the convolution
kernel(set by ModelPsfMatchConfig.kernelSize) is much bigger than: psfSize//2.
Thus the sizes of Psf-model matching kernels are typically smaller than their
image-matching counterparts.  If the size of the kernel is too small, the
convolved stars will look "boxy"; if the kernel is too large, the kernel
solution will be "noisy".  This is a trade-off that needs careful attention
for a given dataset.

The primary use case for this Task is in matching an Exposure to a
constant-across-the-sky Psf model for the purposes of image coaddition. It is
important to note that in the code, the "template" Psf is the Psf that the
science image gets matched to.  In this sense the order of template and
science image are reversed, compared to ImagePsfMatchTask, which operates on
the template image.

.. _lsst.ip.diffim.ModelPsfMatchTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.ModelPsfMatchTask

.. _lsst.ip.diffim.ModelPsfMatchTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.ModelPsfMatchTask

.. _lsst.ip.diffim.ModelPsfMatchTask-config:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.ModelPsfMatchTask

.. _lsst.ip.diffim.ModelPsfMatchTask-debug:

Debugging
=========

The ``pipetask`` command line interface supports a flag ``--debug`` to import
ref:  supports a flag ``--debug`` to import debug.py from your PYTHONPATH;
see :ref:`lsstDebug` for more about debug.py files.
The available variables in ModelsPsfMatchTask include:

display : `bool`
    Enable debug display output.
maskTransparency : `float`
    Transparency of mask planes in the output display.
displaySpatialCells : `bool`
    Show spatial cells before the fit.
