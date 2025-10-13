.. lsst-task-topic:: lsst.ip.diffim.SimplifiedSubtractTask

######################
SimplifiedSubtractTask
######################

.. _lsst.ip.diffim.SimplifiedSubtractTask-description:

Description
===========

This task implements a streamlined image differencing pipeline designed to
operate without requiring an input source catalog, while optionally allowing
the reuse of an existing point spread function (PSF) matching kernel. The
pipeline supports two primary execution methods, each aligned with distinct
scientific objectives:

*Existing Kernel Method*—Utilizes a previously computed PSF matching kernel to
reproduce the results of a prior image differencing run.

*Derived Kernel Method*—Computes the PSF matching kernel internally by performing
a streamlined source detection and measurement process, enabling image
differencing with minimal required inputs.


Simplified Image Differencing Pipeline with the Existing Kernel Method
----------------------------------------------------------------------

Use a local version of the AP pipeline configuration to run the pipeline by reusing a precomputed PSF matching kernel. The sample pipeline file, ``simpleDiffim_existingKernel.yaml``, shown below, demonstrates the recommended configuration settings for this method.

.. code-block:: yaml

    description: A simplified single-visit difference image pipeline that reuses an existing PSF matching kernel to reproduce the results of a previous image differencing run.
    instrument: lsst.obs.lsst.LsstCam
    tasks:
      rewarpTemplate:
        class: lsst.ip.diffim.getTemplate.GetTemplateTask
        config:
          connections.bbox: preliminary_visit_image.bbox
          connections.wcs: preliminary_visit_image.wcs
          connections.coaddExposures: template_coadd
          connections.template: template_detector
      simplifiedSubtractImages:
        class: lsst.ip.diffim.SimplifiedSubtractTask
        config:
          connections.template: template_detector
          connections.science: preliminary_visit_image
          connections.difference: difference_image_predetection
          connections.matchedTemplate: template_matched
          connections.inputPsfMatchingKernel: difference_kernel
          useExistingKernel: True

      detectAndMeasureDiaSource:
        class: lsst.ip.diffim.detectAndMeasure.DetectAndMeasureTask
        config:
          connections.science: preliminary_visit_image
          connections.matchedTemplate: template_matched
          connections.difference:  difference_image_predetection
          connections.outputSchema: dia_source_schema
          connections.diaSources: dia_source_unfiltered
          connections.subtractedMeasuredExposure: difference_image
          connections.maskedStreaks: goodSeeingDiff_streaks
          doSkySources: True
          doCalculateResidualMetics: False


Simplified Image Differencing Pipeline with the Derived Kernel Method
---------------------------------------------------------------------

Alternatively, the pipeline can compute the PSF matching kernel internally by performing source detection and measurement. This method requires the minimal set of inputs necessary to execute image differencing. The sample pipeline tailored for this mode, ``simpleDiffim_derivedKernel.yaml``, shown below, includes the recommended configuration settings for the Derived Kernel Method.

.. code-block:: yaml

    description: A simplified single-visit difference image pipeline to produce image differencing results with the minimal reqired inputs.
    instrument: lsst.obs.lsst.LsstCam
    tasks:
      rewarpTemplate:
        class: lsst.ip.diffim.getTemplate.GetTemplateTask
        config:
          connections.bbox: preliminary_visit_image.bbox
          connections.wcs: preliminary_visit_image.wcs
          connections.coaddExposures: template_coadd
          connections.template: template_detector
      simplifiedSubtractImages:
        class: lsst.ip.diffim.SimplifiedSubtractTask
        config:
          connections.template: template_detector
          connections.science: preliminary_visit_image
          connections.difference: difference_image_predetection
          connections.matchedTemplate: template_matched
          connections.kernelSources: difference_kernel_sources
          connections.psfMatchingKernel: difference_kernel
          useExistingKernel: False

      detectAndMeasureDiaSource:
        class: lsst.ip.diffim.detectAndMeasure.DetectAndMeasureTask
        config:
          connections.science: preliminary_visit_image
          connections.matchedTemplate: template_matched
          connections.difference:  difference_image_predetection
          connections.kernelSources: difference_kernel_sources
          connections.outputSchema: dia_source_schema
          connections.diaSources: dia_source_unfiltered
          connections.subtractedMeasuredExposure: difference_image
          connections.maskedStreaks: goodSeeingDiff_streaks
          doSkySources: True
          doCalculateResidualMetics: True

Execute the Simplified Image Differencing Pipeline Using ``pipetask run``
-------------------------------------------------------------------------

To execute the pipline using the ``pipetask run`` method, implement the following code:

.. code-block:: python

    pipetask run -b /repo/main -i /input/collections -o /your/output/collection -d "instrument='LSSTCam' and exposure=2025050100367 and detector=30 and skymap='lsst_cells_v1'" -p /path/to/your/simpleDiffim.yaml


.. _lsst.ip.diffim.SimplifiedSubtractTask-api:

Python API summary
==================

.. lsst-task-api-summary:: lsst.ip.diffim.SimplifiedSubtractTask

.. _lsst.ip.diffim.SimplifiedSubtractTask-subtasks:

Retargetable subtasks
=====================

.. lsst-task-config-subtasks:: lsst.ip.diffim.SimplifiedSubtractTask

.. _lsst.ip.diffim.SimplifiedSubtractTask-configs:

Configuration fields
====================

.. lsst-task-config-fields:: lsst.ip.diffim.SimplifiedSubtractTask
