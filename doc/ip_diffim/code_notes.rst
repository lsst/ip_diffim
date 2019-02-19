###########################################
ip_diffim package usage and technical notes
###########################################

This page is a collection of usage and code related notes about the
image differencing implmentation. We do not summarise the Alard-Lupton
(AL) [AL_1998]_ and Zackay, Ofek, Gal-Yam (ZOGY) [ZOGY2016]_ papers
themselves here. 

ImageDifferenceTask supports two types of templates that are
subtracted from calexp science exposures: coadds and calexps.  This is
selected by the ``getTemplate`` configurable field. Only the
``getTemplate`` subtask is a properly retargetable top level
``pexConfig.ConfigurableField``.

There is a sequence of pre-subtraction and post-subtraction processing
steps around the actual subtraction of the images. Including the
subtraction operation itself these steps are controlled by the top
level ``do<ACTION>`` configuration options. These top-level
configuration options are summarised in the `subtasks flowchart
<https://github.com/lsst-dm/diffimTests/tree/master/figure_subtasks>`_
. Some of these top level configuration options are passed on to
invoked subtasks, too.

Usage modes
-----------

- template as coadd
  The ``--id`` option specifies a science exposure.
  The ``--templateId`` specifies a skymap.
- template as calexp

- example run with coadd as template
- example run with calexp as template
- tbd example run with dcr coadd as template

  
Using the ``GetCalexpAsTemplate`` subtask, we can select two calexps
for subtraction.  The ``--id`` option specifies a science exposure.
The ``--templateId`` specifies a new *visit* id only. All other fields
are discarded at the ``--templateId`` option; the template
calexp specification inherits all other field values from the ``--id``
option.  The following example subtract visit 410915 ccd 25 (template)
from 411055 ccd 25 (science exposure):
  
.. code-block:: text
	imageDifference.py repo/calibimgs --id visit=411055
		ccd=25 --templateId visit=410915
		
			
General implementation notes
----------------------------

If the two images have the same origin and extent (same WCS) the
subtraction is performed pixel by pixel. Otherwise, the _template_
exposure is warped and resampled into the WCS frame of the science
exposure.

The AL kernel fitting is entirely implemented in C++, an adaptation of
the `hotpants <https://github.com/acbecker/hotpants>`_ C package
by Becker, while the ZOGY algorithm is entirely implemented in
Pyton based on the high level numpy fft functionality.

Beside the top level configuration options, the AL algorithm has its
own set of separate configuration parameters, while the ZOGY algorithm
does not have any algorithm specific parameters.

AL implementation notes
-----------------------


If ``convolveTemplate==False`` the science exposure is convolved and
the image science minus template image difference is multiplied by -1.


The performance of the AL algorithm was studied indetails in the
[Becker_LDM-227]_ report. This study forms the basis of the AL
algorithm default values; the the degree of the polynomial
multiplicator of the Gaussian kernel basis functions (degGauss), the
degree of the polynomial that is fitted to the spatial variaton of the
solution coefficients accross the image (spatialKernelOrder) and the
default detection thresholds (5.5 sigma).

The result of the AL algorithm is a spatially varying set of
coefficients that describes the fitted model convolution kernel to be
applied to the template image before performing per-pixel subtraction
of the images.  Due to noise in the template image, convolving the
template introduces correlation in the noise in the template
image.

The AL algorithm was improved by an additional *afterburner*
decorrelation to remove the noise correlation in the image
difference. The implemented decorrelation method and mathematical
formulae of the decorrelation kernel is summarised and studied in
[Reiss_DMTN-021]_.

ZOGY implementation notes
-------------------------

[ZOGY2016]_ is free from the assumption that the template is noise
free or specially selected by any other means. We simply deal with two
images with different PSFs and noise characteristics (sigma). In the
basic version of the algorithm, the random noise in the pixels are
assumed to be background dominated i.e. uncorrelated between pixels
and independent of the pixel values. Also we assume that the noise has
zero expectation value i.e. the expectation value of the random noise
is already removed.

ZOGY shows that if these assumptions hold, the difference image noise
is also independent and identically distributed over its pixels.


Pre-convolution is not implemented in the ZOGY algorithm. If
``doPreConvolve==True`` the ``S`` the detection likelihood image
(eq. 12 in [ZOGY16]_) and its *corrected variance* ``S_var`` (the
*denominator* of eq. 25 [ZOGY16]_) are calculated following the
*simple* correction steps presented in the paper Section 3.3.

With the variance planes, exactly the same steps are repeated than on
the data planes, only the subtraction step is replaced by addition.

The nan valueas are removed from the science and template images
before Fourier transformations and replaced by the image mean
values. The mask plane UNMASKEDNAN is set for pixels where any of the
two inputs or the difference result is *nan*.

References
----------

.. [AL_1998] Alard, C.; Lupton, Robert H. A Method for Optimal Image
              Subtraction

.. [Reiss_DMTN-021] Reiss J. David, Lupton, Robert H. DMTN-021:
		    Implementation of Image Difference Decorrelation
	      
.. [ZOGY2016] Zackay B., Ofek E. O., Gal-Yam A.,
	      Proper Image Subtraction—Optimal Transient Detection,
	      Photometry, and Hypothesis Testing, 2016, ApJ, 830, 27

.. [Becker_LDM-227] Becker A. et al. LDM-227 Report on Late Winter2013
		    Production: Image Differencing
