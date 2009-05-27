#include <lsst/ip/diffim/ImageSubtract.h>

namespace exceptions = lsst::pex::exceptions; 
namespace logging    = lsst::pex::logging; 
namespace image      = lsst::afw::image;
namespace math       = lsst::afw::math;
namespace detection  = lsst::afw::detection;
namespace diffim     = lsst::ip::diffim;

//
// Constructors
//
template <typename ImageT, typename VarT>
diffim::PsfMatchingFunctor<ImageT, VarT>::PsfMatchingFunctor(
        lsst::afw::math::KernelList<lsst::afw::math::Kernel> const& basisList
    ) :
    _basisList(basisList),
    _background(0.0),
    _backgroundError(0.0),
    _kernel(boost::shared_ptr<lsst::afw::math::Kernel>()),
    _kernelError(boost::shared_ptr<lsst::afw::math::Kernel>())
{}

//
// Public Member Functions
//

template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctor<ImageT, VarT>::reset() {
    //FOR SOME REASON THE KERNEL RESET DOES NOT WORK AND SEG FAULTS
    //this->_background      = 0.;
    //this->_backgroundError = 0.;
    //this->_kernel.reset();
    //this->_kernelError.reset();
}

/** 
 * @brief Create PSF matching kernel
 * @param imageToConvolve image to apply kernel to
 * @param imageToNotConvolve image whose PSF you want to match to
 * @param varianceEstimate estimate of the variance per pixel
 * @param policy read only policy object
 */
template <typename ImageT, typename VarT>
void diffim::PsfMatchingFunctor<ImageT, VarT>::apply(
        Image const& imageToConvolve, 
        Image const& imageToNotConvolve, 
        Variance const& varianceEstimate, 
        lsst::pex::policy::Policy  const& policy
){
    // Make sure you do not overwrite anyone else's kernels
    reset();

    int const nKernelParams = _basisList.size();
    int const nBackgroundParams = 1;
    int const nParams = nKernelParams + nBackgroundParams;

    boost::timer t;
    t.restart();
    
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nParams, nParams);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(nParams);
    
    ImageList convolvedImageList(nKernelParams);
    ImageList::iterator imageIter = convolvedImageList.begin();
    ImageList::iterator const & imageEnd = convolvedImageList.end();
    KernelList::const_iterator basisIter = _basisList.begin();
    KernelList::const_iterator const & basisEnd = _basisList.end();
    
    int const kernelCols = policy.getInt("kernelCols");
    int const kernelRows = policy.getInt("kernelRows");

    image::BBox bbox(image::PointI((*basisIter)->getCtrX(), (*basisIter)->getCtrY()),
            (*imageIter)->getWidth()  - (*basisIter)->getWidth() + 1,
            (*imageIter)->getHeight()  - (*basisIter)->getHeight()  + 1);


    // Create C_ij in the formalism of Alard & Lupton */
    for (; basisIter != basisEnd; ++basisIter, ++imageIter) {
        /*
         * NOTE : we could also *precompute* the entire template image convolved with these functions
         *        and save them somewhere to avoid this step each time.  however, our paradigm is to
         *        compute whatever is needed on the fly.  hence this step here.
         */

        // Ignore buffers around edge of convolved images :
        //
        // If the kernel has width 5, it has center pixel 2.  The first good pixel
        // is the (5-2)=3rd pixel, which is array index 2, and ends up being the
        // index of the central pixel.
        //
        // You also have a buffer of unusable pixels on the other side, numbered
        // width-center-1.  The last good usable pixel is N-width+center+1.
    
        // Example : the kernel is width = 5, center =  2
        //
        //     ---|---|-c-|---|---|
        //          
        //           the image is width = N
        //           convolve this with the kernel, and you get
        //
        //    |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
        //
        //           g = first/last good pixel
        //           x = bad
        // 
        //           the first good pixel is the array index that has the value "center", 2
        //           the last good pixel has array index N-(5-2)+1
        //           eg. if N = 100, you want to use up to index 97
        //               100-3+1 = 98, and the loops use i < 98, meaning the last
        //               index you address is 97.
        Image::Ptr temp(new Image(imageToConvolve.getDimensions()));
        math::convolve(*temp, imageToConvolve, **basisIter, false);
        *imageIter = Image::Ptr(new Image(temp, bbox));
    } 

    Image subToNotConolve(imageToNotConvolve, bbox);
    Variance subVarianceEstimate(varianceEstimate, bbox);

    double time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do basis convolutions : %.2f s", time);

    t.restart();
     
    basisIter = _basisList.begin();
    imageIter = convolvedImageList.begin();

    std::vector<ImgIterator> iteratorList;
    for (; imageIter != imageEnd; ++imageIter) {
        iteratorList.push_back((**imageIter).begin());
    }
    
    ImgIterator toNotConvolveIter= subToNotConvolve.begin();
    VarIterator varianceIter = subVarianceEstimate.begin();

    // Unit test ImageSubtract_1.py should show
    // Image range : 9 9 -> 31 31 : 2804.000000 2798.191162
    logging::TTrace<8>("lsst.ip.diffim.PsfMatchingFunctor.apply",
                       "Image range : %d %d -> %d %d : %f %f",
                       startCol, startRow, endCol, endRow, 
                       0 + *toConvolveIter, 0 + *toNotConvolveIter);

    //loop over all pixels
    int numPixels = bbox.getWidth()*bbox.getHeight();
    for (int pixel=0; pixel < numPixels; ++numPixels) {
        //*****should this next line be toConvolveIter?
        ImageT const ncImage = *toNotConvolveIter;
        double const iVariance = 1.0 / *varianceIter;
            
        // kernel index i
        for (int i = 0; i != nKernelParams; ++i) {
            ImageT const cdImagei = *locatorList[i];
                
            // kernel index j
            for (int j = i; j != nKernelParams; ++j) {
                ImageT const cdImagej  = *locatorList[j];
                M(i, j) += cdImagei*cdImagej*iVariance;
            } 
                
            B(i) += ncImage*cdImagei*iVariance;
                
            // Constant background term; effectively j = kidxj + 1 */
            M(i, nParams-1) += cdImagei*iVariance;
            
            ++locatorList[i];
        } 
            
        // Background term; effectively i = kidxi + 1 
        B(nParams-1) += ncImage*iVariance;
        M(nParams-1, nParams-1) += 1.0*iVariance;
            
        // Step each accessor
        ++toNotConvolveIter;
        ++varianceIter;                    
    } // pixel loop
    
    
    /** @note If we are going to regularize the solution to M, this is the place
     * to do it 
     */
    
    // Fill in rest of M
    for (int i=0; i < nParams; ++i) {
        for (int j=i+1; j < nParams; ++j) {
            M(j, i) = M(i, j);
        }
    }
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to step through pixels : %.2f s", time);
    t.restart();

    //std::cout << "B eigen : " << B << std::endl;

    // To use Cholesky decomposition, the matrix needs to be symmetric (M is, by
    // design) and positive definite.  
    //
    // Eventually put a check in here to make sure its positive definite
    //
    Eigen::VectorXd solution = Eigen::VectorXd::Zero(nParams);;
    if (!( M.ldlt().solve(B, &solution) )) {
        logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                           "Unable to determine kernel via Cholesky LDL^T");
        if (!( M.llt().solve(B, &solution) )) {
            logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                               "Unable to determine kernel via Cholesky LL^T");
            if (!( M.lu().solve(B, &solution) )) {
                logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                                   "Unable to determine kernel via LU");
                // LAST RESORT
                try {
                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
                    Eigen::MatrixXd const& R = solver.eigenvectors();
                    Eigen::VectorXd eigenValues  = solver.eigenvalues().cwise().inverse();
                    
                    solution = R*eigenValues.asDiagonal()*R.transpose()*B;
                } catch (exceptions::Exception& e) {
                    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                                       "Unable to determine kernel via eigen-values");
                    
                    throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution in PsfMatchingFunctor::apply");
                }
            }
        }
    }
    //std::cout << "solution eigen : " << solution << std::endl;
    //return;

    // Estimate of parameter uncertainties comes from the inverse of the
    // covariance matrix (noise spectrum).  
    // N.R. 15.4.8 to 15.4.15
    // 
    // Since this is a linear problem no need to use Fisher matrix
    // N.R. 15.5.8

    // Although I might be able to take advantage of the solution above.
    // Since this now works and is not the rate limiting step, keep as-is for DC3a.

    // Use Cholesky decomposition again.
    // Cholkesy:
    // Cov       =  L L^t
    // Cov^(-1)  = (L L^t)^(-1)
    //           = (L^T)^-1 L^(-1)
    Eigen::MatrixXd covariance = M.transpose() * M;
    Eigen::LLT<Eigen::MatrixXd> llt = covariance.llt();
    Eigen::MatrixXd uncertainty = llt.matrixL().transpose().inverse() * llt.matrixL().inverse();
    
    time = t.elapsed();
    logging::TTrace<5>("lsst.ip.diffim.PsfMatchingFunctor.apply", 
                       "Total compute time to do matrix math : %.2f s", time);
    
    // Translate from Eigen vectors into LSST classes
    std::vector<double> values(kCols*kRows);
    std::vector<double> errValues(kCols*kRows);
    for (int row = 0, idx = 0; row < kernelRows; row++) {
        for (int col = 0; col < kernelCols; col++, idx++) {
            
            // sanity checking
            if (std::isnan(solution(idx))) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel solution (nan)");
            }
            if (std::isnan(uncertainty(idx, idx) )) {
                throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
            }
            if (uncertainty(idx, idx) < 0.0) {
                throw LSST_EXCEPT(exceptions::Exception,
                        str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                                uncertainty(idx, idx)
                        )
                );
            }
            
            values[idx] = solution(idx);
            errValues[idx] = sqrt(uncertainty(idx, idx));
        }
    }
    _kernel.reset(new math::LinearCombinationKernel(_basisList, values));
    _kernelError.reset(new math::LinearCombinationKernel(_basisList, errValues));
    
    // Estimate of Background and Background Error */
    if (std::isnan(uncertainty(nParams-1, nParams-1))) {
        throw LSST_EXCEPT(exceptions::Exception, "Unable to determine kernel uncertainty (nan)");
    }
    if (uncertainty(nParams-1, nParams-1) < 0.0) {
        throw LSST_EXCEPT(exceptions::Exception, 
                str(boost::format("Unable to determine kernel uncertainty, negative variance (%.3e)") % 
                        uncertainty(nParams-1, nParams-1)
                )
        );
    }

    _background = solution(nParams-1);
    _backgroundError = sqrt(uncertainty(nParams-1, nParams-1));
}

