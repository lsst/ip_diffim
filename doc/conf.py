"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.conf.pipelinespkg import *


project = "ip_diffim"
html_theme_options["logotext"] = project
html_title = project
html_short_title = project
doxylink = {}
