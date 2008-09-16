// -*- LSST-C++ -*- // fixed format comment for emacs
/**
 * \file
 *
 * \ingroup diffim
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

#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/exceptions.h>
#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/math/FunctionLibrary.h>
#include <lsst/afw/image/Image.h>
#include <lsst/afw/image/ImageUtils.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/afw/image/MaskedImage.h>
#include <lsst/afw/image/PixelAccessors.h>
#include <lsst/pex/logging/Trace.h> 
#include <lsst/afw/image/Wcs.h> 

#include "lsst/ip/diffim/wcsMatch.h"

/** \brief Remap an Exposure to a new WCS.
 *
 * kernelType is one of:
 * * nearest: nearest neighbor (not yet implemented)
 *   Good noise conservation, bad aliasing issues. Best used for weight maps.
 *   Kernel size must be 3x3
 * * bilinear:
 *   Good for undersampled data fast (not yet implemented)
 *   Kernel size must be 3x3
 * * lanczos: 2-d Lanczos function:
 *   Accurate but slow.
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
 * \throw lsst::pex::exceptions::InvalidParameter error if kernelType is not one of the
 * types listed above.
 * \throw lsst::pex::exceptions::InvalidParameter error if kernelCols != 3 or kernelRows != 3
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
int lsst::ip::diffim::wcsMatch(
    lsst::afw::image::Exposure<ImageT, MaskT> &remapExposure,      ///< remapped exposure
    lsst::afw::image::Exposure<ImageT, MaskT> const &origExposure, ///< original exposure
    std::string const kernelType,   ///< kernel type (see main function docs for more info)
    const int kernelCols,   ///< kernel size - columns
    const int kernelRows    ///< kernel size - rows
    )
{
    int numEdgePixels = 0;

    // Create remapping AnalyticKernel of desired type and size.
    typedef typename lsst::afw::math::AnalyticKernel::KernelFunctionPtr FunctionPtr;
    FunctionPtr akPtrFcn; ///< analytic kernel pointer function

    if (kernelType == "nearest") { 
        // Nearest Neighbor Interpolation. Not yet implemented.
        if (kernelCols != 3 || kernelRows != 3) {
            throw lsst::pex::exceptions::InvalidParameter(
                boost::format("Kernel size must be 3x3 for kernelType %s") % kernelType);
        }
    } else if (kernelType == "bilinear") { 
        // Bi-Linear Interpolation. Not yet implemented.
        if (kernelCols != 3 || kernelRows != 3) {
            throw lsst::pex::exceptions::InvalidParameter(
                boost::format("Kernel size must be 3x3 for kernelType %s") % kernelType);
        }
    } else if (kernelType == "lanczos") { 
        // 2D Lanczos resampling kernel
        int order = (std::min(kernelRows, kernelCols) - 1)/2;
        akPtrFcn = FunctionPtr(new lsst::afw::math::LanczosFunction2<lsst::afw::math::Kernel::PixelT>(order));
    } else {
        throw lsst::pex::exceptions::InvalidParameter(
            boost::format("Invalid kernelType %s") % kernelType);
    }
    if (!akPtrFcn) {
        throw lsst::pex::exceptions::InvalidParameter(
            boost::format("kernelType %s not yet implemented") % kernelType);
    }

    lsst::afw::math::AnalyticKernel remapKernel(*akPtrFcn, kernelCols, kernelRows);
    lsst::pex::logging::Trace("lsst.ip.diffim", 3,
        boost::format("Created analytic kernel of type=%s; cols=%d; rows=%d")
        % kernelType % kernelCols % kernelRows);
    
    // Compute kernel extent; use to prevent applying kernel outside of origExposure
    int colBorder0 = remapKernel.getCtrCol();
    int rowBorder0 = remapKernel.getCtrRow();
    int colBorder1 = remapKernel.getCols() - 1 - colBorder0;
    int rowBorder1 = remapKernel.getRows() - 1 - rowBorder0;

    // Create a blank kernel image of the appropriate size and get the accessor to it.
    lsst::afw::image::Image<lsst::afw::math::Kernel::PixelT> kImage(kernelCols, kernelRows);
    typename vw::ImageView<lsst::afw::math::Kernel::PixelT>::pixel_accessor kAcc = kImage.origin();

    // Get the original MaskedImage and a pixel accessor to it.
    lsst::afw::image::MaskedImage<ImageT, MaskT> origMI = origExposure.getMaskedImage();
    const int origCols = static_cast<int>(origMI.getCols());
    const int origRows = static_cast<int>(origMI.getRows());
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> origMiAcc(origMI);
    lsst::afw::image::Wcs origWcs = origExposure.getWcs();
    lsst::pex::logging::Trace("lsst.ip.diffim", 3,
        boost::format("orig image cols=%d; rows=%d") % origCols % origRows);

    // Get the remapped MaskedImage, a pixel accessor to it and the remapped wcs.
    lsst::afw::image::MaskedImage<ImageT, MaskT> remapMI = remapExposure.getMaskedImage();
    lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> remapRowAcc(remapMI);
    lsst::afw::image::Wcs remapWcs = remapExposure.getWcs();
   
    // Conform mask plane names of remapped MaskedImage to match original
    remapMI.getMask()->conformMaskPlanes(origMI.getMask()->getMaskPlaneDict());
    
    // Make a pixel mask from the EDGE bit, if available (0 if not available)
    const MaskT edgePixelMask = origMI.getMask()->getPlaneBitMask("EDGE");
    lsst::pex::logging::Trace("lsst.ip.diffim", 3, boost::format("edgePixelMask=0x%X") % edgePixelMask);
    
    const int numRemapCols = static_cast<int>(remapMI.getCols());
    const int numRemapRows = static_cast<int>(remapMI.getRows());
    lsst::pex::logging::Trace("lsst.ip.diffim", 3,
        boost::format("remap image cols=%d; rows=%d") % numRemapCols % numRemapRows);

    // The original image accessor points to (0,0) which corresponds to pixel colBorder0, rowBorder0
    // because the accessor points to (0,0) of the kernel rather than the center of the kernel
    int origCol = colBorder0;
    int origRow = rowBorder0;

    // Set each pixel of remapExposure's MaskedImage
    lsst::pex::logging::Trace("lsst.ip.diffim", 4, "Remapping masked image");
    lsst::afw::image::Coord2D remapPosColRow;   
    for (int remapRow = 0; remapRow < numRemapRows; remapRow++, remapRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImageT, MaskT> remapColAcc = remapRowAcc;
        remapPosColRow[1] = lsst::afw::image::indexToPosition(remapRow);
        for (int remapCol = 0; remapCol < numRemapCols; remapCol++, remapColAcc.nextCol()) {
            // compute sky position associated with this pixel of remapped MaskedImage
            remapPosColRow[0] = lsst::afw::image::indexToPosition(remapCol);
            lsst::afw::image::Coord2D raDec = remapWcs.colRowToRaDec(remapPosColRow);            
            
            // compute associated pixel position on original MaskedImage
            lsst::afw::image::Coord2D origPosColRow = origWcs.raDecToColRow(raDec);

            // Compute new corresponding position on original image and break it into integer and fractional
            // parts; the latter is used to compute the remapping kernel.
            std::vector<double> fracOrigPix(2);
            int newOrigCol = lsst::afw::image::positionToIndex(fracOrigPix[0], origPosColRow[0]);
            int newOrigRow = lsst::afw::image::positionToIndex(fracOrigPix[1], origPosColRow[1]);
            
            // Check new position before applying it
            if ((newOrigCol - colBorder0 < 0) || (newOrigCol + colBorder1 >= origCols) 
                || (newOrigRow - rowBorder0 < 0) || (newOrigRow + rowBorder1 >= origRows)) {
                // skip this pixel
                *remapColAcc.image = 0;
                *remapColAcc.variance = 0;
                *remapColAcc.mask = edgePixelMask;
                ++numEdgePixels;
//                lsst::pex::logging::Trace("lsst.ip.diffim", 5, "skipping pixel at remapCol=%d; remapRow=%d",
//                    remapCol, remapRow);
                continue;
            }

            // New original pixel position is usable, advance to it
            origMiAcc.advance(newOrigCol - origCol, newOrigRow - origRow);
            origCol = newOrigCol;
            origRow = newOrigRow;
            lsst::afw::image::Coord2D origColRow;   
            origColRow[0] = origCol;
            origColRow[1] = origRow;

            // Compute new kernel image based on fractional pixel position
            remapKernel.setKernelParameters(fracOrigPix); 
            double kSum;
            remapKernel.computeImage(kImage, kSum, false);
            
            // Determine the intensity multipler due to relative pixel scale and kernel sum
            double multFac = remapWcs.pixArea(remapPosColRow) / (origWcs.pixArea(origColRow) * kSum);
           
            // Apply remapping kernel to original MaskedImage to compute remapped pixel
            lsst::afw::math::apply(remapColAcc, origMiAcc, kAcc, kernelCols, kernelRows);

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
lsst::afw::image::Exposure<ImageT, MaskT> lsst::ip::diffim::wcsMatch(
    int &numEdgePixels, ///< number of pixels that were not computed because they were too close to the edge
                        ///< (or off the edge) of origExposure
    lsst::afw::image::Wcs const &remapWcs,  ///< remapped exposure's WCS
    const int numRemapCols,            ///< remapped exposure size - columns
    const int numRemapRows,            ///< remapped exposure size - rows
    lsst::afw::image::Exposure<ImageT, MaskT> const &origExposure, ///< original exposure 
    std::string const kernelType,   ///< kernel type (see main function docs for more info)
    const int kernelCols,   ///< kernel size - columns
    const int kernelRows    ///< kernel size - rows
    )
{
    lsst::afw::image::MaskedImage<ImageT, MaskT> remapMaskedImage(numRemapCols, numRemapRows);
    lsst::afw::image::Exposure<ImageT, MaskT> remapExposure(remapMaskedImage, remapWcs);

    numEdgePixels = lsst::ip::diffim::wcsMatch(
        remapExposure, origExposure, kernelType, kernelCols, kernelRows); 

    return remapExposure;

} // wcsMatch

/************************************************************************************************************/
//
// Explicit instantiations
//
typedef float imagePixelType;

template
int lsst::ip::diffim::wcsMatch(lsst::afw::image::Exposure<boost::uint16_t, lsst::afw::image::maskPixelType> &remapExposure,
                              lsst::afw::image::Exposure<boost::uint16_t, lsst::afw::image::maskPixelType> const &origExposure,
                              std::string const kernelType,
                              const int kernelCols,
                              const int kernelRows);
template
int lsst::ip::diffim::wcsMatch(lsst::afw::image::Exposure<imagePixelType, lsst::afw::image::maskPixelType> &remapExposure,
                              lsst::afw::image::Exposure<imagePixelType, lsst::afw::image::maskPixelType> const &origExposure,
                              std::string const kernelType,
                              const int kernelCols,
                              const int kernelRows);

template
lsst::afw::image::Exposure<boost::uint16_t, lsst::afw::image::maskPixelType> lsst::ip::diffim::wcsMatch(
    int &numEdgePixels,
    lsst::afw::image::Wcs const &remapWcs,
    const int numRemapCols,       
    const int numRemapRows,       
    lsst::afw::image::Exposure<boost::uint16_t, lsst::afw::image::maskPixelType> const &origExposure,
    std::string const kernelType, 
    const int kernelCols,  
    const int kernelRows);

template
lsst::afw::image::Exposure<imagePixelType, lsst::afw::image::maskPixelType> lsst::ip::diffim::wcsMatch(
    int &numEdgePixels,
    lsst::afw::image::Wcs const &remapWcs,
    const int numRemapCols,       
    const int numRemapRows,       
    lsst::afw::image::Exposure<imagePixelType, lsst::afw::image::maskPixelType> const &origExposure,
    std::string const kernelType, 
    const int kernelCols,  
    const int kernelRows);
//
// Why do we need a double image?
//
template
lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> lsst::ip::diffim::wcsMatch(
    int &numEdgePixels,
    lsst::afw::image::Wcs const &remapWcs,
    const int numRemapCols,       
    const int numRemapRows,       
    lsst::afw::image::Exposure<double, lsst::afw::image::maskPixelType> const &origExposure,
    std::string const kernelType, 
    const int kernelCols,  
    const int kernelRows);
