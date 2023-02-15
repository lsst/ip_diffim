# This file is part of ip_diffim.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["TransiNetInterface"]

import torch

from lsst.meas.transiNet.modelPackages import NNModelPackage
from lsst.afw.image import Exposure


class TransiNetInterface:
    """
    The interface between the LSST AP pipeline and a trained pytorch-based
    end2end TransiNet.
    """

    def __init__(self, model_package_name, package_storage_mode, device='cpu'):
        self.model_package_name = model_package_name
        self.package_storage_mode = package_storage_mode
        self.device = device
        self.init_model()

    def init_model(self):
        """Create and initialize an NN model
        """
        model_package = NNModelPackage(self.model_package_name, self.package_storage_mode)
        self.model = model_package.load(self.device)

        # Put the model in evaluation mode instead of training model.
        self.model.eval()

    def prepare_input(self, template, science):
        """
        Convert inputs from numpy arrays, etc. to a torch.tensor blob.

        Parameters
        ----------
        inputs : `list` [`CutoutInputs`]
            Inputs to be scored.

        Returns
        -------
        blob
            Prepared torch tensor blob to run the model on.
        labels
            Truth labels, concatenated into a single list.
        """

        # Convert each cutout to a torch tensor
        template = torch.from_numpy(template)
        science = torch.from_numpy(science)

        # Stack the components to create a single blob
        blob = torch.stack((template, science), dim=0)
        blob = torch.unsquash(blob, 0)

        return blob

    def infer(self, template, science):
        """Return the score of this cutout.

        Parameters
        ----------
        inputs : `list` [`CutoutInputs`]
            Inputs to be scored.

        Returns
        -------
        difference : `lsst.afw.Exposure`
        """
        blob = self.prepare_input(template, science)
        result = self.model(blob)
        difference = Exposure(result.detach().numpy())

        return difference
