#################################
ip_diffim package technical notes
#################################



In th
The ImageDifferenceTask supports two types of templates that are subtracted from calexp science exposures.
This is set by the configuration options


Usage modes
-----------

- template as coadd
  The ``--id`` option specifies a science exposure.
  The ``--templateId`` specifies a skymap.
- template as calexp

- example run with coadd as template
- example run with calexp as template
- tbd example run with dcr coadd as template

.. code-block:: text
		imageDifference.py repo3/calibimgs --id visit=411055
		ccd=25 --templateId visit=410915
		
  
General implementation notes
----------------------------
It is always the template exposure that is warped and resampled into
the WCS frame of the science exposure.


AL implementation notes
-----------------------


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




If ``convolveTemplate==True`` the science exposure is convolved and the image difference is multiplied by -1.

Pre-convolution is not implemented in the ZOGY algorithm. If
``doPreConvolve==True`` the ``S`` the detection likelihood image
(eq. 12 in [ZOGY16]_) and the its corrected variance ``S_var`` (the
_denominator_ of eq. 25 [ZOGY16]_) are calculated.

With the variance planes, exactly the same steps are repeated than on the data planes, only the subtraction step is replaced by addition.

Test line.

The nan valueas are removed from the science and template images
before Fourier transformations. The mask plane UNMASKEDNAN is set for
pixels where any of the two inputs or the difference result is _nan_.

In the image model of this algorithm, the

.. [ZOGY2016] Zackay B., Ofek E. O., Gal-Yam A.,
	      Proper Image Subtraction—Optimal Transient Detection,
	      Photometry, and Hypothesis Testing, 2016, ApJ, 830, 27
	      
