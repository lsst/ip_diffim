"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.ip.diffim


_g = globals()
_g.update(build_package_configs(
    project_name='ip_diffim',
    version=lsst.ip.diffim.version.__version__))
