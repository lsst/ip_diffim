Example pipelines

*** Image Copy Pipeline

A trivial example of reading a FITS file in and writing it back out.
This pipeline does not use events; instead it informs the Input and Output
stages of the desired file paths by creating input and output policy files
from templates.

To use:
% python runCopyPipeline.py [--help]

*** Image Subtraction Pipeline

Subtract one pair of images specified on the command line.
This pipeline does not use events; instead it informs the Input and Output
stages of the desired file paths by creating input and output policy files
from templates.

To use:
% python runSubtractPipeline.py [--help]

*** Image Many Subtract Pipeline

Subtract a list of pairs of images specified in a text file.
All subtractions are performed on one node, sequentially.
This pipeline uses events and it requires a connection to lsst8.ncsa.uiuc.edu to handle those events.
I plan to remove that requirement once I learn how to configure the pipeline appropriately.

To use:
* Create a file listing the images to subtract.
  For details on the file format run:
  % python feedManySubtractPipeline.p --help
* Start the pipeline:
  % python startManySubtractPipeline.py [--help]
* Wait for the pipeline to start receiving events.
* Feed the pipeline by running the following from another login:
  % python feedDManySubtractPipeline.py [--help]
