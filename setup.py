#!/usr/bin/env python
"""Installer for lsst.

To use:
- Install the prerequisites listed in doc/index.html
- Then do the usual LSST thing:

python setup.py install --lsst-home (or --lsst-devel)

To do:
- Build a recursive system that searches for setup.py files
  and runs them, creating the empty structure and __init__.py files
  as needed. It must have support for building in the right order.
"""
import os
import sys
import glob
from distutils.core import setup
from numarray.numarrayext import NumarrayExtension

PkgBase = "lsst.imageproc"
PyDir = ""


# list all packages and subpackages here
packages = [
	"lsst.imageproc",
	"lsst.imageproc.WCS",
]


# get setuputil
currSysPath = sys.path
#sys.path = [os.path.join(PyDir, "apps", "support")] + list(currSysPath)
import lsst.support.setuputil as setuputil
sys.path = currSysPath

# process sys.argv to handle --lsst-home, etc.
setuputil.procArgv()


print "packages=", packages
setup(
	name = PkgBase,
	description = "LSST Image Processing Pipeline",
	package_dir = {PkgBase: PyDir},
	packages = packages,
)
