// -*- lsst-c++ -*-
/**
 * @file BasisSets.cc
 *
 * @brief Implementation of image subtraction functions declared in BasisSets.h
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <cmath> 
#include <limits>

#include "boost/timer.hpp" 

#include "lsst/pex/exceptions/Exception.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/ip/diffim/BasisSets.h"

namespace pexExcept  = lsst::pex::exceptions; 
namespace afwImage   = lsst::afw::image;
namespace afwMath    = lsst::afw::math;

namespace lsst { 
namespace ip { 
namespace diffim {

/** 
 * @brief Generate a basis set of delta function Kernels.
 *
 * Generates a vector of Kernels sized nCols * nRows, where each Kernel has
 * a unique pixel set to value 1.0 with the other pixels valued 0.0.  This
 * is the "delta function" basis set.
 * 
 * @return Vector of orthonormal delta function Kernels.
 *
 * @throw lsst::pex::exceptions::DomainError if nRows or nCols not positive
 *
 * @ingroup ip_diffim
 */
lsst::afw::math::KernelList
generateDeltaFunctionBasisSet(
    unsigned int width,                 ///< number of columns in the set
    unsigned int height                 ///< number of rows in the set
    ) {
    if ((width < 1) || (height < 1)) {
        throw LSST_EXCEPT(pexExcept::Exception, "nRows and nCols must be positive");
    }
    int const signedWidth = static_cast<int>(width);
    int const signedHeight = static_cast<int>(height);
    afwMath::KernelList kernelBasisList;
    for (int row = 0; row < signedHeight; ++row) {
        for (int col = 0; col < signedWidth; ++col) {
            boost::shared_ptr<afwMath::Kernel> 
                kernelPtr(new afwMath::DeltaFunctionKernel(width, height, afwImage::PointI(col, row)));
            kernelBasisList.push_back(kernelPtr);
        }
    }
    return kernelBasisList;
}

/** 
 * @brief Generate an Alard-Lupton basis set of Kernels.
 *
 * @note Should consider implementing as SeparableKernels for additional speed,
 * but this will make the normalization a bit more complicated
 * 
 * @return Vector of Alard-Lupton Kernels.
 *
 * @ingroup ip_diffim
 */
lsst::afw::math::KernelList
generateAlardLuptonBasisSet(
    unsigned int halfWidth,                ///< size is 2*N + 1
    unsigned int nGauss,                   ///< number of gaussians
    std::vector<double> const &sigGauss,   ///< width of the gaussians
    std::vector<int>    const &degGauss    ///< local spatial variation of gaussians
    ) {
    typedef afwMath::Kernel::Pixel Pixel;
    typedef afwImage::Image<Pixel> Image;

    if (halfWidth < 1) {
        throw LSST_EXCEPT(pexExcept::Exception, "halfWidth must be positive");
    }
    if (nGauss != sigGauss.size()) {
        throw LSST_EXCEPT(pexExcept::Exception, "sigGauss does not have enough entries");
    }
    if (nGauss != degGauss.size()) {
        throw LSST_EXCEPT(pexExcept::Exception, "degGauss does not have enough entries");
    }
    int fullWidth = 2 * halfWidth + 1;
    Image image(fullWidth, fullWidth);
    
    afwMath::KernelList kernelBasisList;
    for (unsigned int i = 0; i < nGauss; i++) {
        /* 
           sigma = FWHM / ( 2 * sqrt(2 * ln(2)) )
        */
        double sig        = sigGauss[i];
        unsigned int deg  = degGauss[i];

        afwMath::GaussianFunction2<Pixel> gaussian(sig, sig);
        afwMath::AnalyticKernel kernel(fullWidth, fullWidth, gaussian);
        afwMath::PolynomialFunction2<Pixel> polynomial(deg);

        for (unsigned int j = 0, n = 0; j <= deg; j++) {
            for (unsigned int k = 0; k <= (deg - j); k++, n++) {
                /* for 0th order term, skip polynomial */
                (void)kernel.computeImage(image, true);
                if (n == 0) {
                    boost::shared_ptr<afwMath::Kernel> 
                        kernelPtr(new afwMath::FixedKernel(image));
                    kernelBasisList.push_back(kernelPtr);
                    continue;
                }
                
                /* gaussian to be modified by this term in the polynomial */
                polynomial.setParameter(n, 1.0);
                (void)kernel.computeImage(image, true);
                for (int y = 0, v = -halfWidth; y < image.getHeight(); y++, v++) {
                    int u = -halfWidth;
                    for (Image::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y); 
                         ptr != end; ++ptr.x(), u++) {
                        /* Evaluate from -1 to 1 */
                        *ptr  = *ptr * polynomial(u/static_cast<double>(halfWidth), 
                                                  v/static_cast<double>(halfWidth));
                    }
                }
                boost::shared_ptr<afwMath::Kernel> 
                    kernelPtr(new afwMath::FixedKernel(image));
                kernelBasisList.push_back(kernelPtr);
                polynomial.setParameter(n, 0.0);
            }
        }
    }
    return renormalizeKernelList(kernelBasisList);
}

/** 
 * @brief Generate regularization matrix for delta function kernels
 */
boost::shared_ptr<Eigen::MatrixXd>
generateFiniteDifferenceRegularization(
    unsigned int width,
    unsigned int height,
    unsigned int order,
    BoundStyle boundaryStyle, // 0 = unwrapped, 1 = wrapped, 2 = order-tapered ('order' is highest used)
    DiffStyle differenceStyle // 0 = forward, 1 = central
    ) {
    
    boost::shared_ptr<Eigen::MatrixXd> bMat = details::generateFdrBMatrix(width, height, order,
                                                                          boundaryStyle, differenceStyle);

    boost::shared_ptr<Eigen::MatrixXd> hMat (new Eigen::MatrixXd(bMat->transpose() * (*bMat)));
    return hMat;

}

namespace details {
/** 
 * @brief Generate regularization matrix for delta function kernels
 */
boost::shared_ptr<Eigen::MatrixXd>
generateFdrBMatrix(
    unsigned int width,
    unsigned int height,
    unsigned int order,
    BoundStyle boundaryStyle, // 0 = unwrapped, 1 = wrapped, 2 = order-tapered ('order' is highest used)
    DiffStyle differenceStyle // 0 = forward, 1 = central
    ) {

    if (order > 2) throw LSST_EXCEPT(pexExcept::Exception, "Only orders 0..2 allowed");

    if ( !(boundaryStyle == UNWRAPPED || boundaryStyle == WRAPPED || boundaryStyle == TAPERED) ) { 
        throw LSST_EXCEPT(pexExcept::Exception, "Boundary styles UNWRAPPED, WRAPPED, TAPERED defined");
    }
    if ( !(differenceStyle == FORWARD_DIFFERENCE || differenceStyle == CENTRAL_DIFFERENCE )) {
        throw LSST_EXCEPT(pexExcept::Exception,
                          "Only FORWARD_DIFFERENCE, and CENTRAL_DIFFERENCE difference types defined");
    }

    /* what works, and what doesn't */
    // == good job == 
    // - order 0, wrapped, forward
    // - order 1, wrapped or tapered, central or forward
    // - order 2, wrapped or tapered, forward
    // == bad job (usually diagongal stripes) ==
    // - all others


    /* 
       Instead of Taylor expanding the forward difference approximation of
       derivatives (N.R. section 18.5) lets just hard code in the expansion of
       the 1st through 3rd derivatives, which will try and enforce smoothness of
       0 through 2nd derivatives.

       A property of the basic "finite difference regularization" is that their
       rows (column multipliers) sum to 0.

       Another consideration is to use *multiple* finite difference operators as
       a constraint.

    */


    // ===========================================================================
    // Get the coeffs for the finite differencing
    // note: The coeffs are stored 2D although they are essentially 1D entities.
    //       The 2d design was chosen to allow cross-terms to be included,
    //         though they are not yet implemented.
    //
    std::vector<std::vector<std::vector<float> > > 
        coeffs(3, std::vector<std::vector<float> >(5, std::vector<float>(5, 0)));
    unsigned int xCen = 0,  yCen = 0;  // center of reqested order coeffs
    unsigned int xCen1 = 0, yCen1 = 0; // center of order 1 coeffs
    unsigned int xCen2 = 0, yCen2 = 0; // center of order 2 coeffs
    unsigned int xSize = 0, ySize = 0;

    // forward difference coefficients
    if (differenceStyle == FORWARD_DIFFERENCE) {
        
        yCen  = xCen  = 0;
        xCen1 = yCen1 = 0;
        xCen2 = yCen2 = 0;

        xSize = ySize = order + 2;

        // default forward difference suggested in NR chap 18
        // 0th order
        coeffs[0][0][0] = -2; 
        coeffs[0][0][1] =  1; 

        coeffs[0][1][0] =  1;  
        coeffs[0][1][1] =  0;

        // 1st 2
        coeffs[1][0][0] = -2; 
        coeffs[1][0][1] =  2;  
        coeffs[1][0][2] = -1; 

        coeffs[1][1][0] =  2;  
        coeffs[1][1][1] =  0;  
        coeffs[1][1][2] =  0; 

        coeffs[1][2][0] = -1; 
        coeffs[1][2][1] =  0;  
        coeffs[1][2][2] =  0; 

        // 2nd 2
        coeffs[2][0][0] = -2; 
        coeffs[2][0][1] =  3;  
        coeffs[2][0][2] = -3; 
        coeffs[2][0][3] =  1; 

        coeffs[2][1][0] =  3;  
        coeffs[2][1][1] =  0;  
        coeffs[2][1][2] =  0; 
        coeffs[2][1][3] =  0; 

        coeffs[2][2][0] = -3; 
        coeffs[2][2][1] =  0;  
        coeffs[2][2][2] =  0; 
        coeffs[2][2][3] =  0; 

        coeffs[2][3][0] =  1;  
        coeffs[2][3][1] =  0;  
        coeffs[2][3][2] =  0; 
        coeffs[2][3][3] =  0; 

    }

    // central difference coefficients
    if (differenceStyle == CENTRAL_DIFFERENCE) {

        // this is asymmetric and produces diagonal banding in the kernel
        // from: http://www.holoborodko.com/pavel/?page_id=239
        if (order == 0) { 
            yCen = xCen = 1;
            xSize = ySize = 3;
        }
        coeffs[0][0][0] =  0; 
        coeffs[0][0][1] = -1; 
        coeffs[0][0][2] =  0; 

        coeffs[0][1][0] = -1; 
        coeffs[0][1][1] =  0; 
        coeffs[0][1][2] =  1; 

        coeffs[0][2][0] =  0; 
        coeffs[0][2][1] =  1; 
        coeffs[0][2][2] =  0; 


        // this works well and is largely the same as order=1 forward-diff.
        // from: http://www.holoborodko.com/pavel/?page_id=239
        if (order == 1) { 
            yCen = xCen = 1;
            xSize = ySize = 3;
        }
        yCen1 = xCen1 = 1;
        coeffs[1][0][0] =  0; 
        coeffs[1][0][1] =  1; 
        coeffs[1][0][2] =  0;  

        coeffs[1][1][0] =  1; 
        coeffs[1][1][1] = -4; 
        coeffs[1][1][2] =  1; 

        coeffs[1][2][0] =  0; 
        coeffs[1][2][1] =  1; 
        coeffs[1][2][2] =  0;  


        // asymmetric and produces diagonal banding in the kernel
        // from http://www.holoborodko.com/pavel/?page_id=239
        if (order == 2) { 
            yCen = xCen = 2;
            xSize = ySize = 5;
        }
        yCen2 = xCen2 = 2;
        coeffs[2][0][0] =  0;
        coeffs[2][0][1] =  0; 
        coeffs[2][0][2] = -1; 
        coeffs[2][0][3] =  0; 
        coeffs[2][0][4] =  0; 

        coeffs[2][1][0] =  0;
        coeffs[2][1][1] =  0; 
        coeffs[2][1][2] =  2; 
        coeffs[2][1][3] =  0; 
        coeffs[2][1][4] =  0; 

        coeffs[2][2][0] = -1;
        coeffs[2][2][1] =  2; 
        coeffs[2][2][2] =  0; 
        coeffs[2][2][3] = -2; 
        coeffs[2][2][4] =  1; 

        coeffs[2][3][0] =  0;
        coeffs[2][3][1] =  0; 
        coeffs[2][3][2] = -2; 
        coeffs[2][3][3] =  0; 
        coeffs[2][3][4] =  0; 

        coeffs[2][4][0] =  0;
        coeffs[2][4][1] =  0; 
        coeffs[2][4][2] =  1; 
        coeffs[2][4][3] =  0; 
        coeffs[2][4][4] =  0; 
    }


    /* Note we have to add 1 extra (empty) term here because of the differential
     * background fitting */
    boost::shared_ptr<Eigen::MatrixXd> bMat(new Eigen::MatrixXd(Eigen::MatrixXd::Zero(width*height + 1,
                                                                                      width*height + 1)));

    /* Forward difference approximation */
    for (unsigned int i = 0; i < width*height; i++) {

        unsigned int const x0 = i % width;  // the x coord in the kernel image
        unsigned int const y0 = i / width;  // the y coord in the kernel image
        
        unsigned int xEdgeDistance = (x0 > (width - x0 - 1))  ? width - x0 - 1  : x0;
        unsigned int yEdgeDistance = (y0 > (height - y0 - 1)) ? height - y0 - 1 : y0;
        unsigned int edgeDistance = (xEdgeDistance < yEdgeDistance) ? xEdgeDistance : yEdgeDistance;

        for (unsigned int dx = 0; dx < xSize; dx++) {
            for (unsigned int dy = 0; dy < ySize; dy++) {

                // determine where to put this coeff

                // handle the boundary condition
                // note: adding width and height in the sum prevents negatives
                unsigned int x = 0;
                unsigned int y = 0; 
                double thisCoeff = 0;

                // no-wrapping at edges
                if (boundaryStyle == UNWRAPPED) {
                    x = x0 + dx - xCen;
                    y = y0 + dy - yCen;
                    if ((y < 0) || (y > height - 1) || (x < 0) || (x > width - 1)) { continue; }
                    thisCoeff = coeffs[order][dx][dy];

                    // wrapping at edges
                } else if (boundaryStyle == WRAPPED) {
                    x = (width  + x0 + dx - xCen) % width;
                    y = (height + y0 + dy - yCen) % height;
                    thisCoeff = coeffs[order][dx][dy];

                    // order tapering to the edge (just clone wrapping for now)
                    // - use the lowest order possible
                } else if (boundaryStyle == TAPERED) {

                    // edge rows and columns ... set to constant
                    if (edgeDistance == 0) {
                        x = x0;
                        y = y0;
                        thisCoeff = 1;
                    }
                    // in one from edge, use 1st order
                    else if (edgeDistance == 1 && order > 0) {
                        x = (width  + x0 + dx - xCen1) % width;
                        y = (height + y0 + dy - yCen1) % height;
                        if ((dx < 3) && (dy < 3)) { thisCoeff = coeffs[1][dx][dy]; } 
                    }
                    // in two from edge, use 2st order if order > 1
                    else if (edgeDistance == 2 && order > 1){
                        x = (width  + x0 + dx - xCen2) % width;
                        y = (height + y0 + dy - yCen2) % height;
                        if ((dx < 5) && (dy < 5)) { thisCoeff = coeffs[2][dx][dy]; } 
                    } 
                    // if we're somewhere in the middle
                    else if (edgeDistance > order) {
                        x = (width  + x0 + dx - xCen) % width;
                        y = (height + y0 + dy - yCen) % height;
                        thisCoeff = coeffs[order][dx][dy];
                    }

                } 

                (*bMat)(i, y*width + x) = thisCoeff;
                
            }

        }

    }

    return bMat;
}
    
} // end namespace details
    

    
/** 
 * @brief Rescale an input set of kernels 
 *
 * @return Vector of renormalized kernels
 *
 * @ingroup ip_diffim
 */
lsst::afw::math::KernelList
renormalizeKernelList(
    lsst::afw::math::KernelList const &kernelListIn ///< Input list to be renormalized
    ) { 
    typedef afwMath::Kernel::Pixel Pixel;
    typedef afwImage::Image<Pixel> Image;
    double kSum;

    /* 
       
    We want all the bases except for the first to sum to 0.0.  This allows
    us to achieve kernel flux conservation (Ksum) across the image since all
    the power will be in the first term, which will not vary spatially.
    
    K(x,y) = Ksum * B_0 + Sum_i : a(x,y) * B_i
    
    To do this, normalize all Kernels to sum = 1. and subtract B_0 from all
    subsequent kenrels.  
    
    To get an idea of the relative contribution of each of these basis
    functions later on down the line, lets also normalize them such that 
    
    Sum(B_i)  == 0.0   *and*
    B_i * B_i == 1.0
    
    For completeness 
    
    Sum(B_0)  == 1.0
    B_0 * B_0 != 1.0
    
    */
    afwMath::KernelList kernelListOut;
    if (kernelListIn.size() == 0) {
        return kernelListOut;
    }

    Image image0(kernelListIn[0]->getDimensions());
    for (unsigned int i = 0; i < kernelListIn.size(); i++) {
        if (i == 0) {
            /* Make sure that it is normalized to kSum 1. */
            (void)kernelListIn[i]->computeImage(image0, true);
            boost::shared_ptr<afwMath::Kernel> 
                kernelPtr(new afwMath::FixedKernel(image0));
            kernelListOut.push_back(kernelPtr);

            continue;
        }

        /* Don't normalize here */
        Image image(kernelListIn[i]->getDimensions());
        (void)kernelListIn[i]->computeImage(image, false);
        
        /* Check the kernel sum; if its close to zero don't do anything */
        kSum = 0.0;
        for (int y = 0; y < image.getHeight(); y++) {
            for (Image::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y); 
                 ptr != end; ++ptr.x()) {
                kSum += *ptr;
            }
        }

        if (fabs(kSum) > std::numeric_limits<double>::epsilon()) {
            image /= kSum;
            image -= image0;
        }


        /* Finally, rescale such that the inner product is 1 */
        kSum = 0.0;
        for (int y = 0; y < image.getHeight(); y++) {
            for (Image::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y); 
                 ptr != end; ++ptr.x()) {
                kSum += *ptr * *ptr;
            }
        }
        image /= std::sqrt(kSum);

        boost::shared_ptr<afwMath::Kernel> 
            kernelPtr(new afwMath::FixedKernel(image));
        kernelListOut.push_back(kernelPtr);
    }
    return kernelListOut;
}

}}} // end of namespace lsst::ip::diffim
