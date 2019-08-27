// -*- lsst-c++ -*-
/**
 * @file BasisLists.cc
 *
 * @brief Implementation of image subtraction functions declared in BasisLists.h
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <cmath>
#include <limits>

#include "boost/timer.hpp"

#include "lsst/pex/exceptions/Exception.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/geom.h"
#include "lsst/log/Log.h"
#include "lsst/ip/diffim/BasisLists.h"

namespace pexExcept  = lsst::pex::exceptions;
namespace geom       = lsst::geom;
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
    makeDeltaFunctionBasisList(
        int width,                 ///< number of columns in the set
        int height                 ///< number of rows in the set
        ) {
        if ((width < 1) || (height < 1)) {
            throw LSST_EXCEPT(pexExcept::Exception, "nRows and nCols must be positive");
        }
        const int signedWidth = static_cast<int>(width);
        const int signedHeight = static_cast<int>(height);
        afwMath::KernelList kernelBasisList;
        for (int row = 0; row < signedHeight; ++row) {
            for (int col = 0; col < signedWidth; ++col) {
                std::shared_ptr<afwMath::Kernel>
                    kernelPtr(new afwMath::DeltaFunctionKernel(width, height, geom::Point2I(col,row)));
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
    makeAlardLuptonBasisList(
        int halfWidth,                ///< size is 2*N + 1
        int nGauss,                   ///< number of gaussians
        std::vector<double> const &sigGauss,   ///< width of the gaussians
        std::vector<int>    const &degGauss    ///< local spatial variation of gaussians
        ) {
        typedef afwMath::Kernel::Pixel Pixel;
        typedef afwImage::Image<Pixel> Image;

        if (halfWidth < 1) {
            throw LSST_EXCEPT(pexExcept::Exception, "halfWidth must be positive");
        }
        if (nGauss != static_cast<int>(sigGauss.size())) {
            throw LSST_EXCEPT(pexExcept::Exception, "sigGauss does not have enough entries");
        }
        if (nGauss != static_cast<int>(degGauss.size())) {
            throw LSST_EXCEPT(pexExcept::Exception, "degGauss does not have enough entries");
        }
        int fullWidth = 2 * halfWidth + 1;
        Image image(geom::Extent2I(fullWidth, fullWidth));

        afwMath::KernelList kernelBasisList;
        for (int i = 0; i < nGauss; i++) {
            /*
               sigma = FWHM / ( 2 * sqrt(2 * ln(2)) )
            */
            double sig  = sigGauss[i];
            int deg     = degGauss[i];

            LOGL_DEBUG("TRACE1.ip.diffim.BasisLists.makeAlardLuptonBasisList",
                       "Gaussian %d : sigma %.2f degree %d", i, sig, deg);

            afwMath::GaussianFunction2<Pixel> gaussian(sig, sig);
            afwMath::AnalyticKernel kernel(fullWidth, fullWidth, gaussian);
            afwMath::PolynomialFunction2<Pixel> polynomial(deg);

            for (int j = 0, n = 0; j <= deg; j++) {
                for (int k = 0; k <= (deg - j); k++, n++) {
                    /* for 0th order term, skip polynomial */
                    (void)kernel.computeImage(image, true);
                    if (n == 0) {
                        std::shared_ptr<afwMath::Kernel>
                            kernelPtr(new afwMath::FixedKernel(image));
                        kernelBasisList.push_back(kernelPtr);
                        continue;
                    }

                    /* gaussian to be modified by this term in the polynomial */
                    polynomial.setParameter(n, 1.);
                    (void)kernel.computeImage(image, true);
                    for (int y = 0, v = -halfWidth; y < image.getHeight(); y++, v++) {
                        int u = -halfWidth;
                        for (Image::xy_locator ptr = image.xy_at(0, y),
                                 end = image.xy_at(image.getWidth(), y);
                             ptr != end; ++ptr.x(), u++) {
                            /* Evaluate from -1 to 1 */
                            *ptr  = *ptr * polynomial(u/static_cast<double>(halfWidth),
                                                      v/static_cast<double>(halfWidth));
                        }
                    }
                    std::shared_ptr<afwMath::Kernel>
                        kernelPtr(new afwMath::FixedKernel(image));
                    kernelBasisList.push_back(kernelPtr);
                    polynomial.setParameter(n, 0.);
                }
            }
        }
        return renormalizeKernelList(kernelBasisList);
    }


    Eigen::MatrixXd makeRegularizationMatrix(
        lsst::daf::base::PropertySet ps
        ) {

        /* NOTES
         *
         * The 3-point first derivative central difference (Laplacian) yields
         * diagonal stripes and should not be used.

         coeffs[0][1] = -1. / 2.;
         coeffs[1][0] = -1. / 2.;
         coeffs[1][2] =  1. / 2.;
         coeffs[2][1] =  1. / 2.;

         * The 5-point second derivative central difference (Laplacian) looks great.
         * For smaller lambdas you need a higher border penalty.  In general, if you
         * decrease the regularization strength you should increase the border
         * penalty to avoid noise in the border pixels.  This has been used to value
         * 10 and things still look OK.

         * The 9-point second derivative central difference (Laplacian with off
         * diagonal terms) also looks great.  Seems to have slightly higher-valued
         * border pixels so make boundary larger if using this.  E.g. 1.5.

         * The forward finite difference, first derivative term works good
         *
         * The forward finite difference, second derivative term is slightly banded
         * in LLC
         *
         * The forward finite difference, third derivative term is highly banded in
         * LLC and should not be used.
         *
         */

        std::string regularizationType = ps.getAsString("regularizationType");
        int width   = ps.getAsInt("kernelSize");
        int height  = ps.getAsInt("kernelSize");
        float borderPenalty  = ps.getAsDouble("regularizationBorderPenalty");
        bool fitForBackground = ps.getAsBool("fitForBackground");

        Eigen::MatrixXd bMat;
        if (regularizationType == "centralDifference") {
            int stencil = ps.getAsInt("centralRegularizationStencil");
            bMat = makeCentralDifferenceMatrix(width, height, stencil, borderPenalty, fitForBackground);
        }
        else if (regularizationType == "forwardDifference") {
            std::vector<int> orders = ps.getAsIntArray("forwardRegularizationOrders");
            bMat = makeForwardDifferenceMatrix(width, height, orders, borderPenalty, fitForBackground);
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "regularizationType not recognized");
        }

        Eigen::MatrixXd hMat = bMat.transpose() * bMat;
        return hMat;
    }

   /**
    * @brief Generate regularization matrix for delta function kernels
    */
    Eigen::MatrixXd makeCentralDifferenceMatrix(
        int width,
        int height,
        int stencil,
        float borderPenalty,
        bool fitForBackground
        ) {

        /* 5- or 9-point stencil to approximate the Laplacian; i.e. this is a second
         * order central finite difference.
         *
         * 5-point approximation ignores off-diagonal elements, and is
         *
         *  f(x-1,y) + f(x+1,y) + f(x,y-1) + f(x,y+1) - 4 f(x,y)
         *
         *   0   1   0
         *   1  -4   1
         *   0   1   0
         *
         *
         * 9-point approximation includes diagonals and is
         *
         *   1   4   1
         *   4  -20  4  / 6
         *   1   4   1
         *
         */

        if (borderPenalty < 0)
            throw LSST_EXCEPT(pexExcept::Exception, "Only border penalty of >= 0 allowed");

        std::vector<std::vector<float> >
            coeffs(3, std::vector<float>(3,0));

        if (stencil == 5) {
            /* http://www.physics.arizona.edu/~restrepo/475B/Notes/source/node51.html
             *
             * This is equivalent to a second order central difference summed along
             * the x-axis, and then summed along the y-axis.
             *
             * http://www.holoborodko.com/pavel/?page_id=239
             *
             */
            coeffs[0][1] =  1.;
            coeffs[1][0] =  1.;
            coeffs[1][1] = -4.;
            coeffs[1][2] =  1.;
            coeffs[2][1] =  1.;
        }
        else if (stencil == 9) {
            /* http://www.physics.arizona.edu/~restrepo/475B/Notes/source/node52.html */
            coeffs[0][0] =   1. / 6.;
            coeffs[0][1] =   4. / 6.;
            coeffs[0][2] =   1. / 6.;
            coeffs[1][0] =   4. / 6.;
            coeffs[1][1] = -20. / 6.;
            coeffs[1][2] =   4. / 6.;
            coeffs[2][0] =   1. / 6.;
            coeffs[2][1] =   4. / 6.;
            coeffs[2][2] =   1. / 6.;
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "Only 5- or 9-point Laplacian stencils allowed");
        }

        int nBgTerms = fitForBackground ? 1 : 0;
        Eigen::MatrixXd bMat = Eigen::MatrixXd::Zero(width * height + nBgTerms, width * height + nBgTerms);

        for (int i = 0; i < width*height; i++) {
            int const x0    = i % width;       // the x coord in the kernel image
            int const y0    = i / width;       // the y coord in the kernel image
            int const distX = width - x0 - 1;  // distance from edge of image
            int const distY = height - y0 - 1; // distance from edge of image

            if ( (x0 > 0) && (y0 > 0) && (distX > 0) && (distY > 0) ) {
                for (int dx = -1; dx < 2; dx += 1) {
                    for (int dy = -1; dy < 2; dy += 1) {
                        bMat(i, i + dx + dy * width) += coeffs[dx+1][dy+1];
                    }
                }
            }
            else {
                bMat(i, i) = borderPenalty;
            }
        }

        if (fitForBackground) {
            /* Last row / col should have no regularization since its the background term */
            if (bMat.col(width*height).sum() != 0.) {
                throw LSST_EXCEPT(pexExcept::Exception, "Error 1 in regularization matrix");
            }
            if (bMat.row(width*height).sum() != 0.) {
                throw LSST_EXCEPT(pexExcept::Exception, "Error 2 in regularization matrix");
            }
        }

        return bMat;
    }

   /**
    * @brief Generate regularization matrix for delta function kernels
    */
    Eigen::MatrixXd makeForwardDifferenceMatrix(
        int width,
        int height,
        std::vector<int> const& orders,
        float borderPenalty,
        bool fitForBackground
        ) {

        /*
           Instead of Taylor expanding the forward difference approximation of
           derivatives (N.R. section 18.5) lets just hard code in the expansion of
           the 1st through 3rd derivatives, which will try and enforce smoothness of
           0 through 2nd derivatives.

           A property of the basic "finite difference regularization" is that their
           rows (column multipliers) sum to 0.

           We taper the order of the constraint as you get close to the boundary.
           The actual perimeter values are left unconstrained but we might want to
           consider giving these the same penalties as the other border pixels,
           which is +1 (since the value gets squared).

        */

        if (borderPenalty < 0)
            throw LSST_EXCEPT(pexExcept::Exception, "Only border penalty of >= 0 allowed");

        std::vector<std::vector<float> >
            coeffs(4, std::vector<float>(4,0));

        // penalize border?  this is along the diagonal and gets squared, so sign does not matter
        coeffs[0][0] = borderPenalty;

        // penalize zeroth derivative
        coeffs[1][0] = -1.;
        coeffs[1][1] = +1.;

        // penalize first derivative
        coeffs[2][0] = -1.;
        coeffs[2][1] = +2.;
        coeffs[2][2] = -1.;

        // penalize second derivative
        coeffs[3][0] = -1.;
        coeffs[3][1] = +3.;
        coeffs[3][2] = -3.;
        coeffs[3][3] = +1.;

        int nBgTerms = fitForBackground ? 1 : 0;
        Eigen::MatrixXd bTot  = Eigen::MatrixXd::Zero(width * height + nBgTerms, width * height + nBgTerms);

        std::vector<int>::const_iterator order;
        for (order = orders.begin(); order != orders.end(); order++) {
            if ((*order < 1) || (*order > 3))
                throw LSST_EXCEPT(pexExcept::Exception, "Only orders 1..3 allowed");

            Eigen::MatrixXd bMatX = Eigen::MatrixXd::Zero(width * height + nBgTerms,
                                                          width * height + nBgTerms);
            Eigen::MatrixXd bMatY = Eigen::MatrixXd::Zero(width * height + nBgTerms,
                                                          width * height + nBgTerms);

            for (int i = 0; i < width*height; i++) {
                int const x0 = i % width;         // the x coord in the kernel image
                int const y0 = i / width;         // the y coord in the kernel image

                int distX       = width - x0 - 1; // distance from edge of image
                int orderToUseX = std::min(distX, *order);
                for (int j = 0; j < orderToUseX+1; j++) {
                    bMatX(i, i + j) = coeffs[orderToUseX][j];
                }

                int distY       = height - y0 - 1; // distance from edge of image
                int orderToUseY = std::min(distY, *order);
                for (int j = 0; j < orderToUseY+1; j++) {
                    bMatY(i, i + j * width) = coeffs[orderToUseY][j];
                }
            }
            bTot += bMatX;
            bTot += bMatY;
        }

        if (fitForBackground) {
            /* Last row / col should have no regularization since its the background term */
            if (bTot.col(width*height).sum() != 0.) {
                throw LSST_EXCEPT(pexExcept::Exception, "Error in regularization matrix");
            }
            if (bTot.row(width*height).sum() != 0.) {
                throw LSST_EXCEPT(pexExcept::Exception, "Error in regularization matrix");
            }
        }

        return bTot;
    }


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
                std::shared_ptr<afwMath::Kernel>
                    kernelPtr(new afwMath::FixedKernel(image0));
                kernelListOut.push_back(kernelPtr);

                continue;
            }

            /* Don't normalize here */
            Image image(kernelListIn[i]->getDimensions());
            (void)kernelListIn[i]->computeImage(image, false);
            /* image.writeFits(str(boost::format("in_k%d.fits") % i)); */

            /* Check the kernel sum; if its close to zero don't do anything */
            kSum = 0.;
            for (int y = 0; y < image.getHeight(); y++) {
                for (Image::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y);
                     ptr != end; ++ptr.x()) {
                    kSum += *ptr;
                }
            }

            /* std::numeric_limits<float>::epsilon() ~ e-7
               std::numeric_limits<double>::epsilon() ~ e-16

               If we end up with 2e-16 kernel sum, this still blows up the kernel values.
               Even though the kernels are double, use the float limits instead
            */
            if (fabs(kSum) > std::numeric_limits<float>::epsilon()) {
                image /= kSum;
                image -= image0;
            }

            /* Finally, rescale such that the inner product is 1 */
            kSum = 0.;
            for (int y = 0; y < image.getHeight(); y++) {
                for (Image::xy_locator ptr = image.xy_at(0, y), end = image.xy_at(image.getWidth(), y);
                     ptr != end; ++ptr.x()) {
                    kSum += *ptr * *ptr;
                }
            }
            image /= std::sqrt(kSum);
            /* image.writeFits(str(boost::format("out_k%d.fits") % i));  */

            std::shared_ptr<afwMath::Kernel>
                kernelPtr(new afwMath::FixedKernel(image));
            kernelListOut.push_back(kernelPtr);
        }
        return kernelListOut;
    }



   /**
    * @brief Generate regularization matrix for delta function kernels
    */
    Eigen::MatrixXd makeFiniteDifferenceRegularizationDeprecated(
        unsigned int width,
        unsigned int height,
        unsigned int order,
        unsigned int boundary_style,   // 0 = unwrapped, 1 = wrapped, 2 = order-tapered
        unsigned int difference_style, // 0 = forward, 1 = central
        bool printB // a debug flag ... remove when done.
        ) {

        if (order > 2) throw LSST_EXCEPT(pexExcept::Exception, "Only orders 0..2 allowed");

        if (boundary_style > 2) {
            throw LSST_EXCEPT(pexExcept::Exception, "Boundary styles 0..2 defined");
        }
        if (difference_style > 1) {
            throw LSST_EXCEPT(pexExcept::Exception, "Only forward(0), and central(1) difference types defined");
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
            coeffs(3, std::vector<std::vector<float> >(5, std::vector<float>(5,0)));
        unsigned int x_cen = 0,  y_cen = 0;  // center of reqested order coeffs
        unsigned int x_cen1 = 0, y_cen1 = 0; // center of order 1 coeffs
        unsigned int x_cen2 = 0, y_cen2 = 0; // center of order 2 coeffs
        unsigned int x_size = 0, y_size = 0;

        // forward difference coefficients
        if (difference_style == 0) {

            y_cen  = x_cen  = 0;
            x_cen1 = y_cen1 = 0;
            x_cen2 = y_cen2 = 0;

            x_size = y_size = order + 2;

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
        if (difference_style == 1) {

            // this is asymmetric and produces diagonal banding in the kernel
            // from: http://www.holoborodko.com/pavel/?page_id=239
            if (order == 0) {
                y_cen = x_cen = 1;
                x_size = y_size = 3;
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
                y_cen = x_cen = 1;
                x_size = y_size = 3;
            }
            y_cen1 = x_cen1 = 1;
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
                y_cen = x_cen = 2;
                x_size = y_size = 5;
            }
            y_cen2 = x_cen2 = 2;
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
        Eigen::MatrixXd bMat = Eigen::MatrixXd::Zero(width*height+1, width*height+1);

        /* Forward difference approximation */
        for (unsigned int i = 0; i < width*height; i++) {

            unsigned int const x0 = i % width;  // the x coord in the kernel image
            unsigned int const y0 = i / width;  // the y coord in the kernel image

            unsigned int x_edge_distance = (x0 > (width - x0 - 1))  ? width - x0 - 1  : x0;
            unsigned int y_edge_distance = (y0 > (height - y0 - 1)) ? height - y0 - 1 : y0;
            unsigned int edge_distance = (x_edge_distance < y_edge_distance) ? x_edge_distance :
                y_edge_distance;

            for (unsigned int dx = 0; dx < x_size; dx++) {
                for (unsigned int dy = 0; dy < y_size; dy++) {

                    // determine where to put this coeff

                    // handle the boundary condition
                    // note: adding width and height in the sum prevents negatives
                    unsigned int x = 0;
                    unsigned int y = 0;
                    double this_coeff = 0;

                    // no-wrapping at edges
                    if (boundary_style == 0) {
                        x = x0 + dx - x_cen;
                        y = y0 + dy - y_cen;
                        if (y > height - 1 || x > width - 1) {
                            continue;
                        }
                        this_coeff = coeffs[order][dx][dy];

                        // wrapping at edges
                    } else if (boundary_style == 1) {
                        x = (width  + x0 + dx - x_cen) % width;
                        y = (height + y0 + dy - y_cen) % height;
                        this_coeff = coeffs[order][dx][dy];

                        // order tapering to the edge (just clone wrapping for now)
                        // - use the lowest order possible
                    } else if (boundary_style == 2) {

                        // edge rows and columns ... set to constant
                        if (edge_distance == 0) {
                            x = x0;
                            y = y0;
                            this_coeff = 1;
                        }
                        // in one from edge, use 1st order
                        else if (edge_distance == 1 && order > 0) {
                            x = (width  + x0 + dx - x_cen1) % width;
                            y = (height + y0 + dy - y_cen1) % height;
                            if ((dx < 3) && (dy < 3)) { this_coeff = coeffs[1][dx][dy]; }
                        }
                        // in two from edge, use 2st order if order > 1
                        else if (edge_distance == 2 && order > 1){
                            x = (width  + x0 + dx - x_cen2) % width;
                            y = (height + y0 + dy - y_cen2) % height;
                            if ((dx < 5) && (dy < 5)) { this_coeff = coeffs[2][dx][dy]; }
                        }
                        // if we're somewhere in the middle
                        else if (edge_distance > order) {
                            x = (width  + x0 + dx - x_cen) % width;
                            y = (height + y0 + dy - y_cen) % height;
                            this_coeff = coeffs[order][dx][dy];
                        }

                    }

                    bMat(i, y*width + x) = this_coeff;

                }

            }

        }

        if (printB)  {
            std::cout << bMat << std::endl;
        }

        Eigen::MatrixXd hMat = bMat.transpose() * bMat;
        return hMat;
    }

}}} // end of namespace lsst::ip::diffim
