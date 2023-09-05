################
Package overview
################

``ip_diffim`` contains pipeline tasks to perform image subtraction and transient source detection on calibrated exposures, template images, and their image differences.
The various implementations of these tasks are described below.

Image subtraction tasks
=======================

The image differencing tasks in this package are implementations of the Alard & Lupton (1998) algorithm (AL) [AL_1998]_, and this can be summarised as applying a kernel determination algorithm in order to match observational properties between two images. 
The method to match images varies according to which task is used, in combination with configuration parameters of such tasks.

 - The pipelines provide the `regular` AL implementation in the `~lsst.ip.diffim.AlardLuptonSubtractTask`, where a kernel is determined to transform the template image; 
 - the `auto-convolution` AL version where the kernel and transformation is applied at either the template or science image, according to its PSF properties --this is when the key parameter `mode` in the configuration  of `~lsst.ip.diffim.AlardLuptonSubtractTask` is set to `auto`; 
 - the `pre-convolution` AL implementation of  `~lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask` that uses an ad-hoc convolution  kernel that is applied to the science image before calculating the optimal transformation kernel.


AlardLuptonSubtractTask
-----------------------

This task is the ``ip_diffim`` implementation of AL algorithm, including further improvements such as de-correlation procedure (introduced in [Reiss_DMTN-021]_).
The task inputs are always a `template` image and a `science` image to subtract, and as well associated source catalogs and image properties information for the science exposure.
It allows to the user to configure whether to apply the de-correlation procedure, and some parameters that affect the catalog used for kernel determination,  background subtraction, etc. 

This task is able as well of performing the mentioned `regular` and `auto-convolution` mode of image subtraction. In the former the `template` image is convolved with a kernel,  and in the alternative the _science_ image is the one being transformed with the convolution kernel. 
Finally, the `auto` setting delegates this decision to the task itself, and so when processing a batch of images some might be processed with `convolveTemplate` and the  rest with the `convolveScience` version. 


AlardLuptonPreconvolveSubtractTask
----------------------------------

The `pre-convolve` task runs a preparation stage before applying the AL algorithm itself, so the determination of the kernel for the `template` image is never a "deconvolution" kernel (negative wing profile, as shown in [Kovacs_DMTN-179]_).

In this stage the `science` image is transformed with a convolution kernel :math:`H`, that is tipically the PSF of this image, becoming a "score image" for detection. 
Once this transformation is applied, standard A&L algorithm derives a kernel for the `template` image so it matches the `science score` image properties, and the subtraction is obtained. \
It is needed to clarify that this result is again a "score" image, and Maximum-Likelihood source detection can be run on it without the standard PSF cross-correlation step.

If one wants to obtain the canonical difference image, then first a deconvolution with :math:`H` is needed. 
Any correlations introduced in the noise properties are later taken into account by means of the standard final de-correlation [Reiss_DMTN-021]_.


Source detection and measurement tasks
======================================

Obtaining difference images is a means to the end of transient and variable source detection. 
For LSST to find the variation in brightness of sources in the sky, it is needed to perform the task of finding the pixels with signal and perform a set of measurements on them.

The main objective is to obtain a catalog of position, shapes and brightness of each identified source, and in addition, sets of flags indicating masked pixels, inaccuracies in the determination of measured parameters, and extra pieces of information about the pixel region related to each source.


DetectAndMeasureTask
--------------------

This task is able to perform maximum-likelihood source identification and determination of photometric parameters from difference images that were obtained with the A&L algorithm implementation in `~lsst.ip.diffim.AlardLuptonSubtractTask`.

In this task, a pre-processing stage is needed where the resulting difference image is cross-correlated with the `science` image PSF profile, in order to obtain a "score image", suitable for maximum-likelihood source detection. 


DetectAndMeasureScoreTask
-------------------------

This task obtains sources and its properties from images resulting from the `~lsst.ip.diffim.AlardLuptonPreconvolveSubtractTask`.

In this case, the score image used comes straight from the subtraction task and no additional image processing is needed.


References
==========

.. [AL_1998] Alard, C.; Lupton, Robert H. A Method for Optimal Image
              Subtraction

.. [Reiss_DMTN-021] Reiss J. David, Lupton, Robert H. DMTN-021:
		    Implementation of Image Difference Decorrelation

.. [Kovacs_DMTN-179] Gabor Kovacs. DMTN-179
        The ZOGY image differencing matching kernel and PSF solutions and their practical implementation issues