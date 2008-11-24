Example pipelines

*** Image Copy Pipeline

A trivial example of reading a FITS file in and writing it back out.
This pipeline does not use events; instead it informs the Input and Output
stages of the desired file paths by creating input and output policy files
from templates.

To use:
% python imageCopyPipeline/runPipeline.py [--help]

*** Image Subtraction Pipeline

Subtract one pair of images specified on the command line.
This pipeline does not use events; instead it informs the Input and Output
stages of the desired file paths by creating input and output policy files
from templates.

To use:
% python imageSubtractPipeline/runPipeline.py [--help]

*** Image Many Subtract Pipeline

Subtract a list of pairs of images specified in a text file.
All subtractions are performed on one node, sequentially.
This pipeline uses events and it requires a connection to lsst8.ncsa.uiuc.edu to handle those events.
I plan to remove that requirement once I learn how to configure the pipeline appropriately.

To use:
* Create a file listing the images to subtract.
  For details on the file format run:
  % python imageManySubtractPipeline/feedPipeline.p --help
* Start the pipeline:
  % python imageManySubtractPipeline/runPipeline.p --help
* Wait for the pipeline to start receiving events.
* Feed the pipeline by running the following from another login:
  % python imageManySubtractPipeline/feedPipeline.p --help

*** Image Subtract and Detect Pipeline

Subtract a list of pairs of images and then run detection on them.
The list of images is supplied in the same way as "imageManySubtractPipeline"
and the pipeline is run in the same way except that you may specify two
policy files: one for image subtraction and one for detection. For details run:
  % python imageSubtractDetectPipeline/runPipeline.p --help
