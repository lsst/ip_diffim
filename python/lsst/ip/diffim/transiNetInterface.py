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

import numpy as np
import torch
import logging
logging.basicConfig(level=logging.DEBUG, format="{name} {levelname}: {message}", style="{")
logger = logging.getLogger("TransiNetInterface")

from lsst.meas.transiNet.modelPackages import NNModelPackage


class TransiNetInterface:
    """
    The interface between the LSST AP pipeline and a trained pytorch-based
    end2end TransiNet.
    """

    def __init__(self, model_package_name, package_storage_mode,
                 # Auto-detect the fastest device available, if not specified
                 device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu"),
                 batch_size=12):
        self.model_package_name = model_package_name
        self.package_storage_mode = package_storage_mode
        self.device = device
        self.batch_size = batch_size
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

    def __forward_pass_crops(self, crops):
        """ Stack the crops into a single blob and pass it through the network.

        Parameters
        ----------
        crops : `list` of `tuple`
            A list of tuples containing the following information:
            (x0, y0, cimg0, cimg1), where (x0, y0) is the position of the
            top-left corner of the crop in the original image, and cimg0
            and cimg1 are the template and science crops, respectively.

        Returns
        -------
        outputs : `torch.Tensor`
            The output of the network.
        """

        # --- create torch blobs out of crops
        blob0 = torch.zeros(
            (len(crops), 1, self.model_input_shape[0], self.model_input_shape[1]), dtype=crops[0][2].dtype)
        blob1 = torch.zeros(
            (len(crops), 1, self.model_input_shape[0], self.model_input_shape[1]), dtype=crops[0][3].dtype)

        for i, (x0, y0, cimg0, cimg1), in enumerate(crops):
            blob0[i, ...] = torch.unsqueeze(cimg0, 0)
            blob1[i, ...] = torch.unsqueeze(cimg1, 0)

        # - Apply the input transform
        blob0 *= self.model_input_scale_factors[0]
        blob1 *= self.model_input_scale_factors[1]

        # - pass the blobs to the network
        i0 = blob0.to(self.device)
        i1 = blob1.to(self.device)

        input = torch.cat((i0, i1), 1)

        output = self.model(input)[0].cpu().detach()

        # - apply the output transform
        output /= self.model_boost_factor

        return output

    def __test_and_distribute_crops(self, crops, output_big, weight_map):
        """ Pass the crops through the network and distribute the results
        back to the output image.

        Parameters
        ----------
        crops : `list` of `tuple`
            A list of tuples containing the following information:
            (x0, y0, cimg0, cimg1), where (x0, y0) is the position of the
            top-left corner of the crop in the original image, and cimg0
            and cimg1 are the template and science crops, respectively.
        output_big : `numpy.ndarray`
            The current output image, with possible results from previous
            crops.
        weight_map : `numpy.ndarray`
            The current weight map, with possible results from previous
            crops.

        Returns
        -------
        output_big : `numpy.ndarray`
            The output image.
        weight_map : `numpy.ndarray`
            The weight map.
        """

        outputs = self.__forward_pass_crops(crops)

        # Distribute the items in the output back to the big-output container
        for i, (x0, y0, cimg0, cimg1), in enumerate(crops):

            output = outputs[i][0].numpy()

            # Apply the window function/mask
            # TODO: this is a hack, we may need to change the window function
            # PREVIOUS CODE: output, mask = applyWindowMask(output)
            mask = np.ones_like(output)

            output_big[y0:self.model_input_shape[1]+y0, x0:self.model_input_shape[0]+x0] += output
            weight_map[y0:self.model_input_shape[1]+y0, x0:self.model_input_shape[0]+x0] += mask

        return output_big, weight_map

    def __test_one_pair(self, img0, img1,
                        batch_size, scanning_stride):
        """ Sweep over the image pair and pass the crops through the network.

        Parameters
        ----------
        img0 : `torch.Tensor`
            The template image.
        img1 : `torch.Tensor`
            The science image.
        batch_size : `int`
            The number of crops to pass through the network at once.
        scanning_stride : `int`
            The number of pixels to move the crop window at each step.

        Returns
        -------
        output_big : `torch.Tensor`
            The output image.
        weight_map : `torch.Tensor`
            The weight map.
        """

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
                if cimg0.shape < tuple(self.model_input_shape[:2]):
                    logger.debug(f'Skipping patch with size: {cimg0.shape} at: {x0}, {y0}')
                    continue

                # Append to the list of crops
                crops.append([x0, y0, cimg0, cimg1])

                # Do a forward pass only if we have gathered enough crops/patches
                if len(crops) == batch_size:
                    output_big, weight_map = self.__test_and_distribute_crops(crops, output_big, weight_map)

                    # --- And free up the memory
                    crops.clear()

        # --- one last round, to take care of the last set, which has had len < batch_size and thus skipped
        if len(crops) > 0:
            output_big, weight_map = self.__test_and_distribute_crops(crops, output_big, weight_map)
            crops.clear()

        # --- undo the effect of sweeping before returning the output
        output_big /= weight_map
        output_big[weight_map == 0] = 0  # we don't want these to be nans, for now

        # --- TODO: Mask-out N/A regions
        # output_big[ref_nan_mask] = 0 #np.nan

        return output_big, weight_map

    def __prepare_input(self, template, science):
        """
        Apply any necessary initial pre-processing to the input images.
        This mainly involves converting the images to torch tensors, but
        does not include cutout-level operations.

        Parameters
        ----------
        template : `lsst.afw.image.Exposure`
            The template image.
        science : `lsst.afw.image.Exposure`
            The science image.

        Returns
        -------
        template : `torch.Tensor`
            The template image as a torch tensor.
        science : `torch.Tensor`
            The science image as a torch tensor.
        """

        # Convert each image to a torch tensor.
        # TODO: see how we want to deal with masks.
        template_ = torch.from_numpy(template.image.array)
        science_ = torch.from_numpy(science.image.array)

        return template_, science_

    def infer(self, template, science):
        """Receive a pair of images and return the difference image
        using TransiNet.

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

        # Convert the images to torch tensors.
        template_, science_ = self.__prepare_input(template, science)

        # Sweep over images, create crops, and pass them through the network.
        scanning_stride = self.model_input_shape[0]//1 - 2*5
        logger.debug(f'Scanning stride: {scanning_stride}')
        output, weight_map = self.__test_one_pair(template_, science_,
                                                  batch_size=self.batch_size,
                                                  scanning_stride=scanning_stride)

        # Convert the result to an Exposure object.
        # TODO: see how we want to deal with masks.
        difference = science.clone()
        difference.image.array = output

        return difference
