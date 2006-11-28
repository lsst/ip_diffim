#!/bin/tcsh 

# Brute force generation of --some-- doco

# Does not handle the documentation for any 'main' codes:
#       Pipeline/WCS/constructWCS/doWCS.py
#       Pipeline/WCS/constructWCS/doMosaicWCS.py
#       Pipeline/WCS/wcsTest/test*.py
# For the moment, just read the top of those files for usage info.

pydoc -w  Pipeline.WCS.constructWCS.WCS  
pydoc -w  Pipeline.WCS.matchCollection.StarSourceMatchCollection  

mkdir -p Doc
mv *.html Doc/
