#################################
ip_diffim package technical notes
#################################

This page is a collection of usage and code related notes about the
image differencing implmentation.

ImageDifferenceTask supports two types of templates that are
subtracted from calexp science exposures: coadds and calexps.  This is
selected by the ``getTemplate`` configurable field. Only the
``getTemplate`` subtask is a properly retargetable top level
``pexConfig.ConfigurableField``. All other subtask execution is controlled by the
top level boolean configuration options.

.. Check which format is necessary for correct object referencing.

This is set by the
configuration options ImageDifferenceTask executes a sequence of
pre-subtraction and post-subtraction steps as well. Including the
subtraction operation itself this is controlled by the top level
``do<ACTION>`` configuration options.


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
The ``--templateId`` specifies a new _visit_ id only. All other fields
are discarded. The template calexp specification inherits all other
field values from the ``--id`` option.  The following example subtract
visit 410915 ccd 25 (template) from 411055 ccd 25 (science exposure):
  
.. code-block:: text
		imageDifference.py repo/calibimgs --id visit=411055
		ccd=25 --templateId visit=410915
		

			
General implementation notes
----------------------------

If the two images have the same origin and extent (same WCS) the
subtraction is performed pixel by pixel. Otherwise, the _template_
exposure is warped and resampled into the WCS frame of the science
exposure.

The task 


AL implementation notes
-----------------------


If ``convolveTemplate==False`` the science exposure is convolved and the image difference is multiplied by -1.


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
(eq. 12 in [ZOGY16]_) and its _corrected_ _variance_ ``S_var`` (the
_denominator_ of eq. 25 [ZOGY16]_) are calculated following the
_simple_ correction steps presented in the paper Section 3.3.

With the variance planes, exactly the same steps are repeated than on the data planes, only the subtraction step is replaced by addition.

Test line.

The nan valueas are removed from the science and template images
before Fourier transformations. The mask plane UNMASKEDNAN is set for
pixels where any of the two inputs or the difference result is _nan_.

In the image model of this algorithm, the

.. [ZOGY2016] Zackay B., Ofek E. O., Gal-Yam A.,
	      Proper Image Subtraction—Optimal Transient Detection,
	      Photometry, and Hypothesis Testing, 2016, ApJ, 830, 27
	      
