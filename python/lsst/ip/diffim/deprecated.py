# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
__all__ = ["deprecate_policy"]

import warnings
import functools
from lsst.pex.policy import Policy


def deprecate_policy(func):
    """Issue a deprecation warning if one of the supplied arguments
    is a Policy, and convert that to a PropertySet."""

    @functools.wraps(func)
    def internal(*args, **kwargs):
        newargs = []
        for i, a in enumerate(args):
            if isinstance(a, Policy):
                warnings.warn(f"pexPolicy in argument {i} is deprecated. Replace with PropertySet"
                              " (Policy support will be removed in v20)",
                              FutureWarning, stacklevel=2)
                a = a.asPropertySet()
            newargs.append(a)
        return func(*newargs, **kwargs)

    return internal
