.. lsst-task-topic:: lsst.ip.diffim.DecorrelateALKernelTask

##########################
DecorrelateALKernelTask
##########################

Pipe-task that removes the neighboring-pixel covariance in an
image difference that are added when the template image is
convolved with the Alard-Lupton PSF matching kernel.

.. _lsst.ip.diffim.DecorrelateALKernelTask-description:

Description
==================

The image differencing pipeline task @link
ip.diffim.psfMatch.PsfMatchTask PSFMatchTask@endlink and @link
ip.diffim.psfMatch.PsfMatchConfigAL PSFMatchConfigAL@endlink uses
the Alard and Lupton (1998) method for matching the PSFs of the
template and science exposures prior to subtraction. The
Alard-Lupton method identifies a matching kernel, which is then
(typically) convolved with the template image to perform PSF
matching. This convolution has the effect of adding covariance
between neighboring pixels in the template image, which is then
added to the image difference by subtraction.

The pixel covariance may be corrected by whitening the noise of
the image difference. This task performs such a decorrelation by
computing a decorrelation kernel (based upon the A&L matching
kernel and variances in the template and science images) and
convolving the image difference with it. This process is described
in detail in [DMTN-021](http://dmtn-021.lsst.io).

This task has no standalone example, however it is applied as a
subtask of pipe.tasks.imageDifference.ImageDifferenceTask.