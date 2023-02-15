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
from lsst.afw.image import ExposureF


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
        """ Create and initialize an NN model and load relevant hyperparameters.
        """
        model_package = NNModelPackage(self.model_package_name, self.package_storage_mode)
        self.model = model_package.load(self.device)

        # Get the hyperparameters
        self.model_input_shape = model_package.get_model_input_shape()
        self.model_input_scale_factors = model_package.get_input_scale_factors()
        self.model_boost_factor = model_package.get_boost_factor()

        # Put the model in evaluation mode
        self.model.eval()

    def forward_pass_crops(self, crops):
        """ Stack the crops into a single blob and pass it through the network.
        """

        # --- create blobs out of crops
        blob0 = np.zeros(
            (len(crops), 1, self.model_input_shape[0], self.model_input_shape[1]), dtype=crops[0][2].dtype)
        blob1 = np.zeros(
            (len(crops), 1, self.model_input_shape[0], self.model_input_shape[1]), dtype=crops[0][2].dtype)
        for i, (x0, y0, cimg0, cimg1), in enumerate(crops):
            blob0[i, ...] = np.expand_dims(cimg0, 0)
            blob1[i, ...] = np.expand_dims(cimg1, 0)

        #- Apply the input transform
        blob0 *= self.model_input_scale_factors[0]
        blob1 *= self.model_input_scale_factors[1]

        # - pass the blobs to the network
        i0 = torch.from_numpy(blob0).float().to(device)
        i1 = torch.from_numpy(blob1).float().to(device)

        input = torch.cat((i0, i1), 1)

        output = self.model(input).cpu().detach().numpy()

        # - apply the output transform
        output /= self.model_boost_factor

        return output

    def test_and_distribute_crops(self, crops, output_big, weight_map):
        
        outputs = forward_pass_crops(crops)
        
        # Distribute the items in the output back to the big-output container
        for i, (x0, y0, cimg0, cimg1), in enumerate(crops):
            
            output = outputs[i][0]
            
            # Apply the window function/mask
            #output, mask = applyWindowMask(output)
            
            output_big[y0:self.model_input_shape[1]+y0, x0:self.model_input_shape[0]+x0] += output
            weight_map[y0:self.model_input_shape[1]+y0, x0:self.model_input_shape[0]+x0] += mask
            
        return output_big, boutput_big, weight_map
    
    def testOnePair(self, img0, img1, ref_nan_mask,
                    batch_size, scanning_stride,
                    results_folder):
        
        # Sweep over the whole image
        output_big = np.zeros_like(img1)
        weight_map = np.zeros_like(output_big)

        # First, generate all the crops and keep them in memory
        crops = []
        for x0 in range(0, img0.shape[1], scanning_stride):
            for y0 in range(0, img0.shape[0], scanning_stride):

                # crop
                cimg0 = img0[y0:self.model_input_shape[1]+y0, x0:self.model_input_shape[0]+x0]
                cimg1 = img1[y0:self.model_input_shape[1]+y0, x0:self.model_input_shape[0]+x0]

                # TODO: decide how to handle near-edge cases
                if cimg0.shape < tuple(self.model_input_shape):
                    print('Skipping patch with size:', cimg0.shape, 'at:', x0, y0)
                    continue

                # Append to the list of crops
                crops.append([x0, y0, cimg0, cimg1])

                # Do a forward pass only if we have gathered enough crops/patches
                if len(crops) == batch_size:
                    output_big, weight_map = self.test_and_distribute_crops(crops, output_big, weight_map)

                    # --- And free up the memory
                    crops.clear()

        # --- one last round, to take care of the last set, which has had len < batch_size and thus skipped
        if len(crops) > 0:
            output_big, weight_map = self.test_and_distribute_crops(crops, output_big, weight_map)
            crops.clear()

        # --- undo the effect of sweeping before returning the output
        output_big /= weight_map
        output_big[weight_map == 0] = 0  # we don't want these to be nans, for now

        # --- Mask-out N/A regions
        # output_big[ref_nan_mask] = 0 #np.nan

        return output_big, weight_map

    # -----------------------------------------------------
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
        # TODO: is getMaskedImage() the right method to use here?
        template = torch.from_numpy(template.getMaskedImage().getImage().getArray())
        science = torch.from_numpy(science.getMaskedImage().getImage().getArray())

        # Stack the components to create a single blob
        blob = torch.stack((template, science), dim=0)
        blob = torch.unsqueeze(blob, dim=0)

        return blob

    def infer(self, template, science):
        """Pass template, science through the model and return the result.

        Parameters
        ----------
        inputs :
            template: `lsst.afw.image.Exposure`
                The template image.
            science: `lsst.afw.image.Exposure`
                The science image.

        Returns
        -------
        difference : `lsst.afw.Exposure`
        """

        # stack the template and science images into a single blob
        blob = self.prepare_input(template, science)

        # Run the model
        result = self.model(blob)[0].detach().numpy().squeeze()

        # Convert the result to an Exposure
        difference = science.clone()
        difference.image.array = result

        return difference
