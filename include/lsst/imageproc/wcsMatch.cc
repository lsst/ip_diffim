// -*- LSST-C++ -*- // fixed format comment for emacs
/**
 * \file
 *
 * \ingroup imageproc
 *
 * \brief Implementation of the templated utility function, wcsMatch, for
 * Astrometric Image Remapping for LSST.  Declared in wcsMatch.h.
 *
 * \author Nicole M. Silvestri, University of Washington
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include <boost/cstdint.hpp> 
#include <boost/format.hpp> 
#include <boost/shared_ptr.hpp>

#include <lsst/mwi/data/DataProperty.h>
#include <lsst/mwi/exceptions.h>
#include <lsst/fw/Exposure.h>
#include <lsst/fw/FunctionLibrary.h>
#include <lsst/fw/Image.h>
#include <lsst/fw/ImageUtils.h>
#include <lsst/fw/Kernel.h>
#include <lsst/fw/KernelFunctions.h>
#include <lsst/fw/MaskedImage.h>
#include <lsst/fw/PixelAccessors.h>
#include <lsst/mwi/utils/Trace.h> 
#include <lsst/fw/WCS.h> 
#include <lsst/imageproc/wcsMatch.h>


/** \brief Remap an Exposure to a new WCS.
 *
 * kernelType is one of:
 * * nearest: nearest neighbor (not yet implemented)
 *   Good noise conservation, bad aliasing issues. Best used for weight maps.
 *   Kernel size must be 3x3
 * * bilinear:
 *   Good for undersampled data fast (not yet implemented)
 *   Kernel size must be 3x3
 * * lanczos: 2-d radial Lanczos function:
 *   Accurate but slow. Not sure if this version or the separable version is best.
 *   f(rad) = sinc(pi rad) sinc(pi rad / n)
 *   with n = (min(kernelRows, kernelCols) - 1)/2
 * * lanczos_separable: 2-d separable Lanczos function:
 *   Accurate but slow. Not sure if this version or the radial version is best.
 *   (x, y) = sinc(pi x') sinc(pi x' / n) sinc(pi y') sinc(pi y' / n)
 *   with n = (min(kernelRows, kernelCols) - 1)/2
 *
 * For pixels in remapExpsure that cannot be computed because their data comes from pixels that are too close
 * to (or off of) the edge of origExposure.
 * * The image and variance are set to 0
 * * The mask bit EDGE is set, if present, else the mask pixel is set to 0
 * * the total number of all such pixels is returned
 *
 * \return the number pixels in remapExpsure that cannot be computed because their data comes from pixels
 * that are too close to (or off of) the edge of origExposure
 *
 * \throw lsst::mwi::exceptions::InvalidParameter error if kernelType is not one of the
 * types listed above.
 * \throw lsst::mwi::exceptions::InvalidParameter error if kernelCols != 3 or kernelRows != 3
 * and kernelType is nearest or bilinear.
 *
 * Algorithm:
 *
 * For each integer pixel position in the remapped Exposure:
 * * The associated sky coordinates are determined using the remapped WCS.
 * * The associated pixel position on origExposure is determined using the original WCS.
 * * A remapping kernel is computed based on the fractional part of the pixel position on origExposure
 * * The remapping kernel is applied to origExposure at the integer portion of the pixel position
 *   to compute the remapped pixel value
 * * The flux-conserving factor is determined from the original and new WCS.
 *   and is applied to the remapped pixel
 *
 * TODO 20071129 Nicole M. Silvestri; By DC3:
 * * Need to synchronize wcsMatch to the UML model robustness/sequence diagrams.
 *   Remove from the Exposure Class in the diagrams.
 *
 * * Should support an additional color-based position correction in the remapping (differential chromatic
 *   refraction). This can be done either object-by-object or pixel-by-pixel.
 *
 * * Need to deal with oversampling and/or weight maps. If done we can use faster kernels than sinc.
 */
template<typename ImageT, typename MaskT> 
int lsst::imageproc::wcsMatch(
    lsst::fw::Exposure<ImageT, MaskT> &remapExposure,      ///< remapped exposure
    lsst::fw::Exposure<ImageT, MaskT> const &origExposure, ///< original exposure
    std::string const kernelType,   ///< kernel type (see main function docs for more info)
    const int kernelCols,   ///< kernel size - columns
    const int kernelRows    ///< kernel size - rows
    )
{
    typedef double KernelType;
    int numEdgePixels = 0;

    // Create remapping AnalyticKernel of desired type and size.
    typedef typename lsst::fw::Kernel<KernelType>::KernelFunctionPtrType funcPtrType;
    funcPtrType akPtrFcn; ///< analytic kernel pointer function

    if (kernelType == "nearest") { 
        // Nearest Neighbor Interpolation. Not yet implemented.
        if (kernelCols != 3 || kernelRows != 3) {
            throw lsst::mwi::exceptions::InvalidParameter(
                boost::format("Kernel size must be 3x3 for kernelType %s") % kernelType);
        }
    } else if (kernelType == "bilinear") { 
        // Bi-Linear Interpolation. Not yet implemented.
        if (kernelCols != 3 || kernelRows != 3) {
            throw lsst::mwi::exceptions::InvalidParameter(
                boost::format("Kernel size must be 3x3 for kernelType %s") % kernelType);
        }
    } else if (kernelType == "lanczos") { 
        // 2D Lanczos resampling kernel: radial version.
        int order = (std::min(kernelRows, kernelCols) - 1)/2; 
        akPtrFcn = funcPtrType(new lsst::fw::function::LanczosFunction2<KernelType>(order));
    } else if (kernelType == "lanczos_separable") { 
        // 2D Lanczos resampling kernel; separable version.
        int order = (std::min(kernelRows, kernelCols) - 1)/2;
        akPtrFcn = funcPtrType(new lsst::fw::function::LanczosSeparableFunction2<KernelType>(order));
    } else {
        throw lsst::mwi::exceptions::InvalidParameter(
            boost::format("Invalid kernelType %s") % kernelType);
    }
    if (!akPtrFcn) {
        throw lsst::mwi::exceptions::InvalidParameter(
            boost::format("kernelType %s not yet implemented") % kernelType);
    }

    lsst::fw::AnalyticKernel<KernelType> remapKernel(akPtrFcn, kernelCols, kernelRows);
    lsst::mwi::utils::Trace("lsst.imageproc", 3,
        boost::format("Created analytic kernel of type=%s; cols=%d; rows=%d")
        % kernelType % kernelCols % kernelRows);
    
    // Compute kernel extent; use to prevent applying kernel outside of origExposure
    int colBorder0 = remapKernel.getCtrCol();
    int rowBorder0 = remapKernel.getCtrRow();
    int colBorder1 = remapKernel.getCols() - 1 - colBorder0;
    int rowBorder1 = remapKernel.getRows() - 1 - rowBorder0;

    // Create a blank kernel image of the appropriate size and get the accessor to it.
    lsst::fw::Image<KernelType> kImage(kernelCols, kernelRows);
    typename vw::ImageView<KernelType>::pixel_accessor kAcc = kImage.origin();

    // Get the original MaskedImage and a pixel accessor to it.
    lsst::fw::MaskedImage<ImageT, MaskT> origMI = origExposure.getMaskedImage();
    const int origCols = static_cast<int>(origMI.getCols());
    const int origRows = static_cast<int>(origMI.getRows());
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> origMiAcc(origMI);
    lsst::fw::WCS origWcs = origExposure.getWcs();
    lsst::mwi::utils::Trace("lsst.imageproc", 3,
        boost::format("orig image cols=%d; rows=%d") % origCols % origRows);

    // Get the remapped MaskedImage, a pixel accessor to it and the remapped wcs.
    lsst::fw::MaskedImage<ImageT, MaskT> remapMI = remapExposure.getMaskedImage();
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> remapRowAcc(remapMI);
    lsst::fw::WCS remapWcs = remapExposure.getWcs();
   
    // Conform mask plane names of remapped MaskedImage to match original
    remapMI.getMask()->conformMaskPlanes(origMI.getMask()->getMaskPlaneDict());
    
    // Make a pixel mask from the EDGE bit, if available (0 if not available)
    const MaskT edgePixelMask = origMI.getMask()->getPlaneBitMask("EDGE");
    lsst::mwi::utils::Trace("lsst.imageproc", 3, boost::format("edgePixelMask=0x%X") % edgePixelMask);
    
    const int numRemapCols = static_cast<int>(remapMI.getCols());
    const int numRemapRows = static_cast<int>(remapMI.getRows());
    lsst::mwi::utils::Trace("lsst.imageproc", 3,
        boost::format("remap image cols=%d; rows=%d") % numRemapCols % numRemapRows);

    // The original image accessor points to (0,0) which corresponds to pixel colBorder0, rowBorder0
    // because the accessor points to (0,0) of the kernel rather than the center of the kernel
    int origCol = colBorder0;
    int origRow = rowBorder0;

    // Set each pixel of remapExposure's MaskedImage
    lsst::mwi::utils::Trace("lsst.imageproc", 4, "Remapping masked image");
    lsst::fw::Coord2D remapPosColRow;   
    for (int remapRow = 0; remapRow < numRemapRows; remapRow++, remapRowAcc.nextRow()) {
        lsst::fw::MaskedPixelAccessor<ImageT, MaskT> remapColAcc = remapRowAcc;
        remapPosColRow[1] = lsst::fw::image::indexToPosition(remapRow);
        for (int remapCol = 0; remapCol < numRemapCols; remapCol++, remapColAcc.nextCol()) {
            // compute sky position associated with this pixel of remapped MaskedImage
            remapPosColRow[0] = lsst::fw::image::indexToPosition(remapCol);
            lsst::fw::Coord2D raDec = remapWcs.colRowToRaDec(remapPosColRow);            
            
            // compute associated pixel position on original MaskedImage
            lsst::fw::Coord2D origPosColRow = origWcs.raDecToColRow(raDec);

            // Compute new corresponding position on original image and break it into integer and fractional
            // parts; the latter is used to compute the remapping kernel.
            std::vector<double> fracOrigPix(2);
            int newOrigCol = lsst::fw::image::positionToIndex(fracOrigPix[0], origPosColRow[0]);
            int newOrigRow = lsst::fw::image::positionToIndex(fracOrigPix[1], origPosColRow[1]);
            
            // Check new position before applying it
            if ((newOrigCol - colBorder0 < 0) || (newOrigCol + colBorder1 >= origCols) 
                || (newOrigRow - rowBorder0 < 0) || (newOrigRow + rowBorder1 >= origRows)) {
                // skip this pixel
                *remapColAcc.image = 0;
                *remapColAcc.variance = 0;
                *remapColAcc.mask = edgePixelMask;
                ++numEdgePixels;
//                lsst::mwi::utils::Trace("lsst.imageproc", 5, "skipping pixel at remapCol=%d; remapRow=%d",
//                    remapCol, remapRow);
                continue;
            }

            // New original pixel position is usable, advance to it
            origMiAcc.advance(newOrigCol - origCol, newOrigRow - origRow);
            origCol = newOrigCol;
            origRow = newOrigRow;
            lsst::fw::Coord2D origColRow;   
            origColRow[0] = origCol;
            origColRow[1] = origRow;

            // Compute new kernel image based on fractional pixel position
            remapKernel.setKernelParameters(fracOrigPix); 
            double kSum;
            remapKernel.computeImage(kImage, kSum, 0.0, 0.0, false);
            
            // Determine the intensity multipler due to relative pixel scale and kernel sum
            double multFac = remapWcs.pixArea(remapPosColRow) / (origWcs.pixArea(origColRow) * kSum);
           
            // Apply remapping kernel to original MaskedImage to compute remapped pixel
            lsst::fw::kernel::apply<ImageT, MaskT, KernelType>(
                remapColAcc, origMiAcc, kAcc, kernelCols, kernelRows);

            // multiply the output from apply function by the computed gain here
            *remapColAcc.image *= static_cast<ImageT>(multFac);
            *remapColAcc.variance *= static_cast<ImageT>(multFac * multFac);

        } // for remapCol loop
    } // for remapRow loop      
    return numEdgePixels;
} // wcsMatch


/** 
 * \brief Remap an Exposure to a new WCS.  
 *
 * This version takes a remapped Exposure's WCS (probably a copy of an existing WCS), the requested size
 * of the remapped MaskedImage, the exposure to be remapped (original Exposure), and the remapping kernel
 * information (kernel type and size).
 *
 * \return the final remapped Exposure
 */
template<typename ImageT, typename MaskT> 
lsst::fw::Exposure<ImageT, MaskT> lsst::imageproc::wcsMatch(
    int &numEdgePixels, ///< number of pixels that were not computed because they were too close to the edge
                        ///< (or off the edge) of origExposure
    lsst::fw::WCS const &remapWcs,  ///< remapped exposure's WCS
    const int numRemapCols,            ///< remapped exposure size - columns
    const int numRemapRows,            ///< remapped exposure size - rows
    lsst::fw::Exposure<ImageT, MaskT> const &origExposure, ///< original exposure 
    std::string const kernelType,   ///< kernel type (see main function docs for more info)
    const int kernelCols,   ///< kernel size - columns
    const int kernelRows    ///< kernel size - rows
    )
{
    lsst::fw::MaskedImage<ImageT, MaskT> remapMaskedImage(numRemapCols, numRemapRows);
    lsst::fw::Exposure<ImageT, MaskT> remapExposure(remapMaskedImage, remapWcs);

    numEdgePixels = lsst::imageproc::wcsMatch(
        remapExposure, origExposure, kernelType, kernelCols, kernelRows); 

    return remapExposure;

} // wcsMatch

// Explicit instantiations

// template void lsst::imageproc::wcsMatch<float, lsst::fw::maskPixelType>;
// template void lsst::imageproc::wcsMatch<double, lsst::fw::maskPixelType>; 
// template void lsst::imageproc::wcsMatch<boost::uint16_t, lsst::fw::maskPixelType>;

