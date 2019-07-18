// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

/**
 * @file ImageSubtract.h
 *
 * @brief Image Subtraction helper functions
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_IMAGESUBTRACT_H
#define LSST_IP_DIFFIM_IMAGESUBTRACT_H

#include "Eigen/Core"

#include "lsst/afw/math.h"
#include "lsst/afw/image.h"

namespace lsst {
namespace ip {
namespace diffim {


    /**
     * @brief Execute fundamental task of convolving template and subtracting it from science image
     *
     * @note This version accepts a MaskedImage for the template
     *
     * @param templateImage  MaskedImage to apply convolutionKernel to
     * @param scienceMaskedImage  MaskedImage from which convolved templateImage is subtracted
     * @param convolutionKernel  Kernel to apply to templateImage
     * @param background  Background scalar or function to subtract after convolution
     * @param invert  Invert the output difference image
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT, typename BackgroundT>
    lsst::afw::image::MaskedImage<PixelT> convolveAndSubtract(
        lsst::afw::image::MaskedImage<PixelT> const& templateImage,
        lsst::afw::image::MaskedImage<PixelT> const& scienceMaskedImage,
        lsst::afw::math::Kernel const& convolutionKernel,
        BackgroundT background,
        bool invert=true
        );

    /**
     * @brief Execute fundamental task of convolving template and subtracting it from science image
     *
     * @note This version accepts an Image for the template, and is thus faster during convolution
     *
     * @param templateImage  Image to apply convolutionKernel to
     * @param scienceMaskedImage  MaskedImage from which convolved templateImage is subtracted
     * @param convolutionKernel  Kernel to apply to templateImage
     * @param background  Background scalar or function to subtract after convolution
     * @param invert  Invert the output difference image
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT, typename BackgroundT>
    lsst::afw::image::MaskedImage<PixelT> convolveAndSubtract(
        lsst::afw::image::Image<PixelT> const& templateImage,
        lsst::afw::image::MaskedImage<PixelT> const& scienceMaskedImage,
        lsst::afw::math::Kernel const& convolutionKernel,
        BackgroundT background,
        bool invert=true
        );

    /**
     * @brief Turns a 2-d Image into a 2-d Eigen Matrix
     *
     * @param img  Image whose pixel values are read into an Eigen::MatrixXd
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    Eigen::MatrixXd imageToEigenMatrix(
        lsst::afw::image::Image<PixelT> const& img
        );

    Eigen::MatrixXi maskToEigenMatrix(
        lsst::afw::image::Mask<lsst::afw::image::MaskPixel> const& mask
        );

}}} // end of namespace lsst::ip::diffim

#endif
