// -*- lsst-c++ -*-
/**
 * @file KernelSolution.cc
 *
 * @brief Implementation of KernelSolution class
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */
#include <iterator>
#include <cmath>
#include <algorithm>
#include <limits>

#include "boost/timer.hpp" 

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/QR"
#include "Eigen/LU"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

#include "lsst/afw/math.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/KernelSolution.h"

#include "lsst/ndarray.h"
#include "lsst/ndarray/eigen.h"

#define DEBUG_MATRIX  0
#define DEBUG_MATRIX2 0

namespace afwDet         = lsst::afw::detection;
namespace afwMath        = lsst::afw::math;
namespace afwGeom        = lsst::afw::geom;
namespace afwImage       = lsst::afw::image;
namespace pexLog         = lsst::pex::logging; 
namespace pexExcept      = lsst::pex::exceptions; 
namespace ndarray        = lsst::ndarray;

namespace lsst { 
namespace ip { 
namespace diffim {
    
    /* Unique identifier for solution */
    int KernelSolution::_SolutionId = 0;

    KernelSolution::KernelSolution(
        boost::shared_ptr<Eigen::MatrixXd> mMat,
        boost::shared_ptr<Eigen::VectorXd> bVec,
        bool fitForBackground
        ) :
        _id(++_SolutionId),
        _mMat(mMat),
        _bVec(bVec),
        _aVec(),
        _solvedBy(NONE),
        _fitForBackground(fitForBackground)
    {};

    KernelSolution::KernelSolution(
        bool fitForBackground
        ) :
        _id(++_SolutionId),
        _mMat(),
        _bVec(),
        _aVec(),
        _solvedBy(NONE),
        _fitForBackground(fitForBackground)
    {};

    KernelSolution::KernelSolution() :
        _id(++_SolutionId),
        _mMat(),
        _bVec(),
        _aVec(),
        _solvedBy(NONE),
        _fitForBackground(true)
    {};

    void KernelSolution::solve() {
        solve(*_mMat, *_bVec);
    }

    double KernelSolution::getConditionNumber(ConditionNumberType conditionType) {
        return getConditionNumber(*_mMat, conditionType);
    }

    double KernelSolution::getConditionNumber(Eigen::MatrixXd mMat, 
                                              ConditionNumberType conditionType) {
        switch (conditionType) {
        case EIGENVALUE: 
            {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(mMat);
            Eigen::VectorXd eValues = eVecValues.eigenvalues();
            double eMax = eValues.maxCoeff();
            double eMin = eValues.minCoeff();
            pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.getConditionNumber", 
                              "EIGENVALUE eMax / eMin = %.3e", eMax / eMin);
            return (eMax / eMin);
            break;
            }
        case SVD: 
            {
            Eigen::VectorXd sValues = mMat.jacobiSvd().singularValues();
            double sMax = sValues.maxCoeff();
            double sMin = sValues.minCoeff();
            pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.getConditionNumber", 
                              "SVD eMax / eMin = %.3e", sMax / sMin);
            return (sMax / sMin);
            break;
            }
        default:
            {
            throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                              "Undefined ConditionNumberType : only EIGENVALUE, SVD allowed.");
            break;
            }
        }
    }

    void KernelSolution::solve(Eigen::MatrixXd mMat,
                               Eigen::VectorXd bVec) {
        
        if (DEBUG_MATRIX) {
            std::cout << "M " << std::endl;
            std::cout << mMat << std::endl;
            std::cout << "B " << std::endl;
            std::cout << bVec << std::endl;
        }

        Eigen::VectorXd aVec = Eigen::VectorXd::Zero(bVec.size());

        boost::timer t;
        t.restart();

        pexLog::TTrace<2>("lsst.ip.diffim.KernelSolution.solve", 
                          "Solving for kernel");
		_solvedBy = LU;
		Eigen::FullPivLU<Eigen::MatrixXd> lu(mMat);
		if (lu.isInvertible()) {
			aVec = lu.solve(bVec);
		} else {
			pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.solve", 
							  "Unable to determine kernel via LU");
			/* LAST RESORT */
			try {
				
				_solvedBy = EIGENVECTOR;
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(mMat);
				Eigen::MatrixXd const& rMat = eVecValues.eigenvectors();
				Eigen::VectorXd eValues = eVecValues.eigenvalues();
				
				for (int i = 0; i != eValues.rows(); ++i) {
					if (eValues(i) != 0.0) {
						eValues(i) = 1.0/eValues(i);
					}
				}
				
				aVec = rMat * eValues.asDiagonal() * rMat.transpose() * bVec;
			} catch (pexExcept::Exception& e) {
				
				_solvedBy = NONE;
				pexLog::TTrace<5>("lsst.ip.diffim.KernelSolution.solve", 
								  "Unable to determine kernel via eigen-values");
				
				throw LSST_EXCEPT(pexExcept::Exception, "Unable to determine kernel solution");
			}
		}

        double time = t.elapsed();
        pexLog::TTrace<3>("lsst.ip.diffim.KernelSolution.solve", 
                          "Compute time for matrix math : %.2f s", time);

        if (DEBUG_MATRIX) {
		  std::cout << "A " << std::endl;
            std::cout << aVec << std::endl;
        }
        
        _aVec = boost::shared_ptr<Eigen::VectorXd>(new Eigen::VectorXd(aVec));
    }

    /*******************************************************************************************************/

    template <typename InputT>
    StaticKernelSolution<InputT>::StaticKernelSolution(
        lsst::afw::math::KernelList const& basisList,
        bool fitForBackground
        ) 
        :
        KernelSolution(fitForBackground),
        _cMat(),
        _iVec(),
        _ivVec(),
        _kernel(),
        _background(0.0),
        _kSum(0.0)
    {
        std::vector<double> kValues(basisList.size());
        _kernel = boost::shared_ptr<afwMath::Kernel>( 
            new afwMath::LinearCombinationKernel(basisList, kValues) 
            );
    };

    template <typename InputT>
    lsst::afw::math::Kernel::Ptr StaticKernelSolution<InputT>::getKernel() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return solution");
        }
        return _kernel;
    }

    template <typename InputT>
    KernelSolution::ImageT::Ptr StaticKernelSolution<InputT>::makeKernelImage() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return image");
        }
        ImageT::Ptr image (
            new ImageT(_kernel->getDimensions())
            );
        (void)_kernel->computeImage(*image, false);              
        return image;
    }

    template <typename InputT>
    double StaticKernelSolution<InputT>::getBackground() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return background");
        }
        return _background;
    }

    template <typename InputT>
    double StaticKernelSolution<InputT>::getKsum() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return ksum");
        }
        return _kSum;
    }

    template <typename InputT>
    std::pair<boost::shared_ptr<lsst::afw::math::Kernel>, double>
    StaticKernelSolution<InputT>::getSolutionPair() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return solution");
        }

        return std::make_pair(_kernel, _background);
    }

    template <typename InputT>
    void StaticKernelSolution<InputT>::build(
        lsst::afw::image::Image<InputT> const &imageToConvolve,
        lsst::afw::image::Image<InputT> const &imageToNotConvolve,
        lsst::afw::image::Image<lsst::afw::image::VariancePixel> const &varianceEstimate
        ) {

        afwMath::Statistics varStats = afwMath::makeStatistics(varianceEstimate, afwMath::MIN);
        if (varStats.getValue(afwMath::MIN) < 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance less than 0.0");
        }
        if (varStats.getValue(afwMath::MIN) == 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance equals 0.0, cannot inverse variance weight");
        }

        lsst::afw::math::KernelList basisList = 
            boost::shared_dynamic_cast<afwMath::LinearCombinationKernel>(_kernel)->getKernelList();
        
        unsigned int const nKernelParameters     = basisList.size();
        unsigned int const nBackgroundParameters = _fitForBackground ? 1 : 0;
        unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;

        std::vector<boost::shared_ptr<afwMath::Kernel> >::const_iterator kiter = basisList.begin();
        
        /* Ignore buffers around edge of convolved images :
         * 
         * If the kernel has width 5, it has center pixel 2.  The first good pixel
         * is the (5-2)=3rd pixel, which is array index 2, and ends up being the
         * index of the central pixel.
         * 
         * You also have a buffer of unusable pixels on the other side, numbered
         * width-center-1.  The last good usable pixel is N-width+center+1.
         * 
         * Example : the kernel is width = 5, center = 2
         * 
         * ---|---|-c-|---|---|
         * 
         * the image is width = N
         * convolve this with the kernel, and you get
         * 
         * |-x-|-x-|-g-|---|---| ... |---|---|-g-|-x-|-x-|
         * 
         * g = first/last good pixel
         * x = bad
         * 
         * the first good pixel is the array index that has the value "center", 2
         * the last good pixel has array index N-(5-2)+1
         * eg. if N = 100, you want to use up to index 97
         * 100-3+1 = 98, and the loops use i < 98, meaning the last
         * index you address is 97.
         */

        /* NOTE - we are accessing particular elements of Eigen arrays using
           these coordinates, therefore they need to be in LOCAL coordinates.
           This was written before ndarray unification.
        */
        afwGeom::Box2I goodBBox = (*kiter)->shrinkBBox(imageToConvolve.getBBox(afwImage::LOCAL));
        unsigned int const startCol = goodBBox.getMinX();
        unsigned int const startRow = goodBBox.getMinY();
        // endCol/Row is one past the index of the last good col/row
        unsigned int endCol = goodBBox.getMaxX() + 1;
        unsigned int endRow = goodBBox.getMaxY() + 1;

        boost::timer t;
        t.restart();
        
        /* Eigen representation of input images; only the pixels that are unconvolved in cimage below */
        Eigen::MatrixXd eigenToConvolve = imageToEigenMatrix(imageToConvolve).block(startRow, 
                                                                                    startCol, 
                                                                                    endRow-startRow, 
                                                                                    endCol-startCol);
        Eigen::MatrixXd eigenToNotConvolve = imageToEigenMatrix(imageToNotConvolve).block(startRow, 
                                                                                          startCol, 
                                                                                          endRow-startRow, 
                                                                                          endCol-startCol);
        Eigen::MatrixXd eigeniVariance = imageToEigenMatrix(varianceEstimate).block(
	    startRow, startCol, endRow-startRow, endCol-startCol
	).array().inverse().matrix();

        /* Resize into 1-D for later usage */
        eigenToConvolve.resize(eigenToConvolve.rows()*eigenToConvolve.cols(), 1);
        eigenToNotConvolve.resize(eigenToNotConvolve.rows()*eigenToNotConvolve.cols(), 1);
        eigeniVariance.resize(eigeniVariance.rows()*eigeniVariance.cols(), 1);
        
        /* Holds image convolved with basis function */
        afwImage::Image<PixelT> cimage(imageToConvolve.getDimensions());
        
        /* Holds eigen representation of image convolved with all basis functions */
        std::vector<boost::shared_ptr<Eigen::MatrixXd> > convolvedEigenList(nKernelParameters);
        
        /* Iterators over convolved image list and basis list */
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiter = 
            convolvedEigenList.begin();
        /* Create C_i in the formalism of Alard & Lupton */
        for (kiter = basisList.begin(); kiter != basisList.end(); ++kiter, ++eiter) {
            afwMath::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */

            boost::shared_ptr<Eigen::MatrixXd> cMat (
                new Eigen::MatrixXd(imageToEigenMatrix(cimage).block(startRow, 
                                                                     startCol, 
                                                                     endRow-startRow, 
                                                                     endCol-startCol))
                );
            cMat->resize(cMat->rows()*cMat->cols(), 1);
            *eiter = cMat;

        } 

        double time = t.elapsed();
        pexLog::TTrace<5>("lsst.ip.diffim.StaticKernelSolution.build", 
                          "Total compute time to do basis convolutions : %.2f s", time);
        t.restart();
        
        /* 
           Load matrix with all values from convolvedEigenList : all images
           (eigeniVariance, convolvedEigenList) must be the same size
        */
        Eigen::MatrixXd cMat(eigenToConvolve.col(0).size(), nParameters);
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterj = 
            convolvedEigenList.begin();
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterE = 
            convolvedEigenList.end();
        for (unsigned int kidxj = 0; eiterj != eiterE; eiterj++, kidxj++) {
            cMat.col(kidxj) = (*eiterj)->col(0);
        }
        /* Treat the last "image" as all 1's to do the background calculation. */
        if (_fitForBackground)
            cMat.col(nParameters-1).fill(1.);

        _cMat.reset(new Eigen::MatrixXd(cMat));
        _ivVec.reset(new Eigen::VectorXd(eigeniVariance.col(0)));
        _iVec.reset(new Eigen::VectorXd(eigenToNotConvolve.col(0)));

        /* Make these outside of solve() so I can check condition number */
        _mMat.reset(new Eigen::MatrixXd((*_cMat).transpose() * ((*_ivVec).asDiagonal() * (*_cMat))));
        _bVec.reset(new Eigen::VectorXd((*_cMat).transpose() * ((*_ivVec).asDiagonal() * (*_iVec))));
    }

    template <typename InputT>
    void StaticKernelSolution<InputT>::solve() {
        pexLog::TTrace<5>("lsst.ip.diffim.StaticKernelSolution.solve", 
                          "mMat is %d x %d; bVec is %d; cMat is %d x %d; vVec is %d; iVec is %d", 
                          (*_mMat).rows(), (*_mMat).cols(), (*_bVec).size(),
                          (*_cMat).rows(), (*_cMat).cols(), (*_ivVec).size(), (*_iVec).size());

        /* If I put this here I can't check for condition number before solving */
        /*
        _mMat.reset(new Eigen::MatrixXd((*_cMat).transpose() * ((*_ivVec).asDiagonal() * (*_cMat))));
        _bVec.reset(new Eigen::VectorXd((*_cMat).transpose() * ((*_ivVec).asDiagonal() * (*_iVec))));
        */

        if (DEBUG_MATRIX) {
            std::cout << "C" << std::endl;
            std::cout << (*_cMat) << std::endl;
            std::cout << "iV" << std::endl;
            std::cout << (*_ivVec) << std::endl;
            std::cout << "I" << std::endl;
            std::cout << (*_iVec) << std::endl;
        }

        try {
            KernelSolution::solve();
        } catch (pexExcept::Exception &e) {
            LSST_EXCEPT_ADD(e, "Unable to solve static kernel matrix");
            throw e;
        }
        /* Turn matrices into _kernel and _background */
        _setKernel();
    }

    template <typename InputT>
    void StaticKernelSolution<InputT>::_setKernel() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot make solution");
        }

        unsigned int const nParameters           = _aVec->size();
        unsigned int const nBackgroundParameters = _fitForBackground ? 1 : 0;
        unsigned int const nKernelParameters     = 
            boost::shared_dynamic_cast<afwMath::LinearCombinationKernel>(_kernel)->getKernelList().size();
        if (nParameters != (nKernelParameters + nBackgroundParameters)) 
            throw LSST_EXCEPT(pexExcept::Exception, "Mismatched sizes in kernel solution");

        /* Fill in the kernel results */
        std::vector<double> kValues(nKernelParameters);
        for (unsigned int idx = 0; idx < nKernelParameters; idx++) {
            if (std::isnan((*_aVec)(idx))) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  str(boost::format("Unable to determine kernel solution %d (nan)") % idx));
            }
            kValues[idx] = (*_aVec)(idx);
        }
        _kernel->setKernelParameters(kValues);

        ImageT::Ptr image (
            new ImageT(_kernel->getDimensions())
            );
        _kSum  = _kernel->computeImage(*image, false);              
        
        if (_fitForBackground) {
            if (std::isnan((*_aVec)(nParameters-1))) {
                throw LSST_EXCEPT(pexExcept::Exception, 
                                  str(boost::format("Unable to determine background solution %d (nan)") % 
                                      (nParameters-1)));
            }
            _background = (*_aVec)(nParameters-1);
        }
    }        


    template <typename InputT>
    void StaticKernelSolution<InputT>::_setKernelUncertainty() {
        throw LSST_EXCEPT(pexExcept::Exception, "Uncertainty calculation not supported");

        /* Estimate of parameter uncertainties comes from the inverse of the
         * covariance matrix (noise spectrum).  
         * N.R. 15.4.8 to 15.4.15
         * 
         * Since this is a linear problem no need to use Fisher matrix
         * N.R. 15.5.8
         *
         * Although I might be able to take advantage of the solution above.
         * Since this now works and is not the rate limiting step, keep as-is for DC3a.
         *
         * Use Cholesky decomposition again.
         * Cholkesy:
         * Cov       =  L L^t
         * Cov^(-1)  = (L L^t)^(-1)
         *           = (L^T)^-1 L^(-1)
         * 
         * Code would be:
         * 
         * Eigen::MatrixXd             cov    = (*_mMat).transpose() * (*_mMat);
         * Eigen::LLT<Eigen::MatrixXd> llt    = cov.llt();
         * Eigen::MatrixXd             error2 = llt.matrixL().transpose().inverse()*llt.matrixL().inverse();
         */
    }

    /*******************************************************************************************************/

    template <typename InputT>
    MaskedKernelSolution<InputT>::MaskedKernelSolution(
        lsst::afw::math::KernelList const& basisList,
        bool fitForBackground
        ) :
        StaticKernelSolution<InputT>(basisList, fitForBackground)
    {};

    template <typename InputT>
    void MaskedKernelSolution<InputT>::build(
        lsst::afw::image::Image<InputT> const &imageToConvolve,
        lsst::afw::image::Image<InputT> const &imageToNotConvolve,
        lsst::afw::image::Image<lsst::afw::image::VariancePixel> const &varianceEstimate,
        lsst::afw::image::Mask<lsst::afw::image::MaskPixel> const &pixelMask
        ) {

        afwMath::Statistics varStats = afwMath::makeStatistics(varianceEstimate, afwMath::MIN);
        if (varStats.getValue(afwMath::MIN) < 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance less than 0.0");
        }
        if (varStats.getValue(afwMath::MIN) == 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance equals 0.0, cannot inverse variance weight");
        }

        /* Full footprint of all input images */
        afwDet::Footprint::Ptr fullFp(new afwDet::Footprint(imageToConvolve.getBBox(afwImage::PARENT)));

        afwMath::KernelList basisList = 
            boost::shared_dynamic_cast<afwMath::LinearCombinationKernel>(this->_kernel)->getKernelList();
        std::vector<boost::shared_ptr<afwMath::Kernel> >::const_iterator kiter = basisList.begin();

        /* Only BAD pixels marked in this mask */
        afwImage::MaskPixel bitMask = 
            (afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("BAD") | 
             afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("SAT") |
             afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("EDGE"));

        /* Create a Footprint that contains all the masked pixels set above */
        afwDet::Threshold threshold = afwDet::Threshold(bitMask, afwDet::Threshold::BITMASK, true);
        afwDet::FootprintSet<InputT, afwImage::MaskPixel> maskFpSet(pixelMask, threshold, true);

        /* And spread it by the kernel half width */
        int growPix = (*kiter)->getCtr().getX();
        afwDet::FootprintSet<InputT, afwImage::MaskPixel> maskedFpSetGrown = 
            afwDet::FootprintSet<InputT, afwImage::MaskPixel>(maskFpSet, growPix, true);
        
        /*
        for (typename afwDet::FootprintSet<InputT, afwImage::MaskPixel>::FootprintList::iterator 
                 ptr = maskedFpSetGrown.getFootprints()->begin(),
                 end = maskedFpSetGrown.getFootprints()->end(); 
             ptr != end; 
             ++ptr) {
            
            afwDet::setMaskFromFootprint(finalMask, 
                                         (**ptr),
                                         afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("BAD"));
        }
        */

        afwImage::Mask<afwImage::MaskPixel> finalMask(pixelMask.getDimensions());
        afwDet::setMaskFromFootprintList(&finalMask, 
                                         *(maskedFpSetGrown.getFootprints()),
                                         afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("BAD"));
        pixelMask.writeFits("pixelmask.fits");
        finalMask.writeFits("finalmask.fits");


        ndarray::Array<int, 1, 1> maskArray = 
            lsst::ndarray::allocate(lsst::ndarray::makeVector(fullFp->getArea()));
        afwDet::flattenArray(*fullFp, finalMask.getArray(), 
                             maskArray, imageToConvolve.getXY0()); /* Need to fake the XY0 */
        ndarray::EigenView<int, 1, 1> maskEigen(maskArray);

        ndarray::Array<InputT, 1, 1> arrayToConvolve = 
            lsst::ndarray::allocate(lsst::ndarray::makeVector(fullFp->getArea()));
        afwDet::flattenArray(*fullFp, imageToConvolve.getArray(), 
                             arrayToConvolve, imageToConvolve.getXY0());
        ndarray::EigenView<InputT, 1, 1> eigenToConvolve0(arrayToConvolve);

        ndarray::Array<InputT, 1, 1> arrayToNotConvolve = 
            lsst::ndarray::allocate(lsst::ndarray::makeVector(fullFp->getArea()));
        afwDet::flattenArray(*fullFp, imageToNotConvolve.getArray(), 
                             arrayToNotConvolve, imageToNotConvolve.getXY0());
        ndarray::EigenView<InputT, 1, 1> eigenToNotConvolve0(arrayToNotConvolve);

        ndarray::Array<afwImage::VariancePixel, 1, 1> arrayVariance = 
            lsst::ndarray::allocate(lsst::ndarray::makeVector(fullFp->getArea()));
        afwDet::flattenArray(*fullFp, varianceEstimate.getArray(), 
                             arrayVariance, varianceEstimate.getXY0());
        ndarray::EigenView<afwImage::VariancePixel, 1, 1> eigenVariance0(arrayVariance);

        int nGood = 0;
        for (int i = 0; i < maskEigen.size(); i++) {
            if (maskEigen(i) == 0.0)
                nGood += 1;
        }

        Eigen::VectorXd eigenToConvolve(nGood);
        Eigen::VectorXd eigenToNotConvolve(nGood);
        Eigen::VectorXd eigenVariance(nGood);
        int nUsed = 0;
        for (int i = 0; i < maskEigen.size(); i++) {
            if (maskEigen(i) == 0.0) {
                eigenToConvolve(nUsed) = eigenToConvolve0(i);
                eigenToNotConvolve(nUsed) = eigenToNotConvolve0(i);
                eigenVariance(nUsed) = eigenVariance0(i);
                nUsed += 1;
            }
        }


        boost::timer t;
        t.restart();

        unsigned int const nKernelParameters     = basisList.size();
        unsigned int const nBackgroundParameters = this->_fitForBackground ? 1 : 0;
        unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;

        /* Holds image convolved with basis function */
        afwImage::Image<InputT> cimage(imageToConvolve.getDimensions());
        
        /* Holds eigen representation of image convolved with all basis functions */
        std::vector<boost::shared_ptr<Eigen::VectorXd> >
            convolvedEigenList(nKernelParameters);
        
        /* Iterators over convolved image list and basis list */
        typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiter = 
            convolvedEigenList.begin();

        /* Create C_i in the formalism of Alard & Lupton */
        for (kiter = basisList.begin(); kiter != basisList.end(); ++kiter, ++eiter) {
            afwMath::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */

            ndarray::Array<InputT, 1, 1> arrayC = 
                lsst::ndarray::allocate(lsst::ndarray::makeVector(fullFp->getArea()));
            afwDet::flattenArray(*fullFp, cimage.getArray(), 
                                 arrayC, cimage.getXY0());
            ndarray::EigenView<InputT, 1, 1> eigenC0(arrayC);

            Eigen::VectorXd eigenC(nGood);
            int nUsed = 0;
            for (int i = 0; i < maskEigen.size(); i++) {
                if (maskEigen(i) == 0.0) {
                    eigenC(nUsed) = eigenC0(i);
                    nUsed += 1;
                }
            }

            boost::shared_ptr<Eigen::VectorXd> 
                eigenCptr (new Eigen::VectorXd(eigenC));
            
            *eiter = eigenCptr;
        }
        double time = t.elapsed();
        pexLog::TTrace<5>("lsst.ip.diffim.StaticKernelSolution.build", 
                          "Total compute time to do basis convolutions : %.2f s", time);
        t.restart();
        
        /* Load matrix with all convolved images */
        Eigen::MatrixXd cMat(eigenToConvolve.size(), nParameters);
        typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterj = 
            convolvedEigenList.begin();
        typename std::vector<boost::shared_ptr<Eigen::VectorXd> >::iterator eiterE = 
            convolvedEigenList.end();
        for (unsigned int kidxj = 0; eiterj != eiterE; eiterj++, kidxj++) {
            //cMat.col(kidxj) = (*eiterj)->template cast<double>();
            cMat.col(kidxj) = (**eiterj);
        }
        /* Treat the last "image" as all 1's to do the background calculation. */
        if (this->_fitForBackground)
            cMat.col(nParameters-1).fill(1.);
        
        this->_cMat.reset(new Eigen::MatrixXd(cMat));
        //this->_ivVec.reset(new Eigen::VectorXd((eigenVariance.template cast<double>()).cwise().inverse()));
        //this->_iVec.reset(new Eigen::VectorXd(eigenToNotConvolve.template cast<double>()));
        this->_ivVec.reset(new Eigen::VectorXd(eigenVariance.array().inverse().matrix()));
        this->_iVec.reset(new Eigen::VectorXd(eigenToNotConvolve));

        /* Make these outside of solve() so I can check condition number */
        this->_mMat.reset(
            new Eigen::MatrixXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_cMat))
            );
        this->_bVec.reset(
            new Eigen::VectorXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_iVec))
            );
    }


    template <typename InputT>
    void MaskedKernelSolution<InputT>::buildOrig(
        lsst::afw::image::Image<InputT> const &imageToConvolve,
        lsst::afw::image::Image<InputT> const &imageToNotConvolve,
        lsst::afw::image::Image<lsst::afw::image::VariancePixel> const &varianceEstimate,
        lsst::afw::image::Mask<lsst::afw::image::MaskPixel> pixelMask
        ) {

        afwMath::Statistics varStats = afwMath::makeStatistics(varianceEstimate, afwMath::MIN);
        if (varStats.getValue(afwMath::MIN) < 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance less than 0.0");
        }
        if (varStats.getValue(afwMath::MIN) == 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance equals 0.0, cannot inverse variance weight");
        }

        lsst::afw::math::KernelList basisList = 
            boost::shared_dynamic_cast<afwMath::LinearCombinationKernel>(this->_kernel)->getKernelList();
        std::vector<boost::shared_ptr<afwMath::Kernel> >::const_iterator kiter = basisList.begin();

        /* Only BAD pixels marked in this mask */
        lsst::afw::image::Mask<lsst::afw::image::MaskPixel> sMask(pixelMask, true);
        afwImage::MaskPixel bitMask = 
            (afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("BAD") | 
             afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("SAT") |
             afwImage::Mask<afwImage::MaskPixel>::getPlaneBitMask("EDGE"));
        sMask &= bitMask;
        /* TBD: Need to figure out a way to spread this; currently done in Python */

        unsigned int const nKernelParameters     = basisList.size();
        unsigned int const nBackgroundParameters = this->_fitForBackground ? 1 : 0;
        unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;

        /* NOTE - we are accessing particular elements of Eigen arrays using
           these coordinates, therefore they need to be in LOCAL coordinates.
           This was written before ndarray unification.
        */
        /* Ignore known EDGE pixels for speed */
        afwGeom::Box2I shrunkLocalBBox = (*kiter)->shrinkBBox(imageToConvolve.getBBox(afwImage::LOCAL));
        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Limits of good pixels after convolution: %d,%d -> %d,%d (local)", 
                          shrunkLocalBBox.getMinX(), shrunkLocalBBox.getMinY(), 
                          shrunkLocalBBox.getMaxX(), shrunkLocalBBox.getMaxY());

        /* Subimages are addressed in raw pixel coordinates */
        unsigned int startCol = shrunkLocalBBox.getMinX();
        unsigned int startRow = shrunkLocalBBox.getMinY();
        unsigned int endCol   = shrunkLocalBBox.getMaxX();
        unsigned int endRow   = shrunkLocalBBox.getMaxY();
        /* Needed for Eigen block slicing operation */
        endCol += 1;
        endRow += 1;

        boost::timer t;
        t.restart();

        /* Eigen representation of input images; only the pixels that are unconvolved in cimage below */
        Eigen::MatrixXi eMask = maskToEigenMatrix(sMask).block(startRow, 
                                                               startCol, 
                                                               endRow-startRow, 
                                                               endCol-startCol);

        Eigen::MatrixXd eigenToConvolve = imageToEigenMatrix(imageToConvolve).block(startRow, 
                                                                                    startCol, 
                                                                                    endRow-startRow, 
                                                                                    endCol-startCol);
        Eigen::MatrixXd eigenToNotConvolve = imageToEigenMatrix(imageToNotConvolve).block(startRow, 
                                                                                          startCol, 
                                                                                          endRow-startRow, 
                                                                                          endCol-startCol);
        Eigen::MatrixXd eigeniVariance = imageToEigenMatrix(varianceEstimate).block(
	    startRow, startCol, endRow-startRow, endCol-startCol
	).array().inverse().matrix();

        /* Resize into 1-D for later usage */
        eMask.resize(eMask.rows()*eMask.cols(), 1);
        eigenToConvolve.resize(eigenToConvolve.rows()*eigenToConvolve.cols(), 1);
        eigenToNotConvolve.resize(eigenToNotConvolve.rows()*eigenToNotConvolve.cols(), 1);
        eigeniVariance.resize(eigeniVariance.rows()*eigeniVariance.cols(), 1);

        /* Do the masking, slowly for now... */
        Eigen::MatrixXd maskedEigenToConvolve(eigenToConvolve.rows(), 1);
        Eigen::MatrixXd maskedEigenToNotConvolve(eigenToNotConvolve.rows(), 1);
        Eigen::MatrixXd maskedEigeniVariance(eigeniVariance.rows(), 1);
        int nGood = 0;
        for (int i = 0; i < eMask.rows(); i++) {
            if (eMask(i, 0) == 0) {
                maskedEigenToConvolve(nGood, 0) = eigenToConvolve(i, 0);
                maskedEigenToNotConvolve(nGood, 0) = eigenToNotConvolve(i, 0);
                maskedEigeniVariance(nGood, 0) = eigeniVariance(i, 0);
                nGood += 1;
            }
        }
        /* Can't resize() since all values are lost; use blocks */
        eigenToConvolve    = maskedEigenToConvolve.block(0, 0, nGood, 1);
        eigenToNotConvolve = maskedEigenToNotConvolve.block(0, 0, nGood, 1);
        eigeniVariance     = maskedEigeniVariance.block(0, 0, nGood, 1);


        /* Holds image convolved with basis function */
        afwImage::Image<InputT> cimage(imageToConvolve.getDimensions());
        
        /* Holds eigen representation of image convolved with all basis functions */
        std::vector<boost::shared_ptr<Eigen::MatrixXd> > convolvedEigenList(nKernelParameters);
        
        /* Iterators over convolved image list and basis list */
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiter = 
            convolvedEigenList.begin();
        /* Create C_i in the formalism of Alard & Lupton */
        for (kiter = basisList.begin(); kiter != basisList.end(); ++kiter, ++eiter) {
            afwMath::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */
            
            Eigen::MatrixXd cMat = imageToEigenMatrix(cimage).block(startRow, 
                                                                    startCol, 
                                                                    endRow-startRow, 
                                                                    endCol-startCol);
            cMat.resize(cMat.rows() * cMat.cols(), 1);

            /* Do masking */
            Eigen::MatrixXd maskedcMat(cMat.rows(), 1);
            int nGood = 0;
            for (int i = 0; i < eMask.rows(); i++) {
                if (eMask(i, 0) == 0) {
                    maskedcMat(nGood, 0) = cMat(i, 0);
                    nGood += 1;
                }
            }
            cMat = maskedcMat.block(0, 0, nGood, 1);
            boost::shared_ptr<Eigen::MatrixXd> cMatPtr (new Eigen::MatrixXd(cMat));
            *eiter = cMatPtr;
        } 

        double time = t.elapsed();
        pexLog::TTrace<5>("lsst.ip.diffim.StaticKernelSolution.build", 
                          "Total compute time to do basis convolutions : %.2f s", time);
        t.restart();
        
        /* 
           Load matrix with all values from convolvedEigenList : all images
           (eigeniVariance, convolvedEigenList) must be the same size
        */
        Eigen::MatrixXd cMat(eigenToConvolve.col(0).size(), nParameters);
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterj = 
            convolvedEigenList.begin();
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterE = 
            convolvedEigenList.end();
        for (unsigned int kidxj = 0; eiterj != eiterE; eiterj++, kidxj++) {
            cMat.col(kidxj) = (*eiterj)->col(0);
        }
        /* Treat the last "image" as all 1's to do the background calculation. */
        if (this->_fitForBackground)
            cMat.col(nParameters-1).fill(1.);

        this->_cMat.reset(new Eigen::MatrixXd(cMat));
        this->_ivVec.reset(new Eigen::VectorXd(eigeniVariance.col(0)));
        this->_iVec.reset(new Eigen::VectorXd(eigenToNotConvolve.col(0)));

        /* Make these outside of solve() so I can check condition number */
        this->_mMat.reset(
            new Eigen::MatrixXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_cMat))
            );
        this->_bVec.reset(
            new Eigen::VectorXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_iVec))
            );
        
    }

    /* NOTE - this was written before the ndarray unification.  I am rewriting
       buildSingleMask to use the ndarray stuff */
    template <typename InputT>
    void MaskedKernelSolution<InputT>::buildSingleMaskOrig(
        lsst::afw::image::Image<InputT> const &imageToConvolve,
        lsst::afw::image::Image<InputT> const &imageToNotConvolve,
        lsst::afw::image::Image<lsst::afw::image::VariancePixel> const &varianceEstimate,
        lsst::afw::geom::Box2I maskBox
        ) {

        afwMath::Statistics varStats = afwMath::makeStatistics(varianceEstimate, afwMath::MIN);
        if (varStats.getValue(afwMath::MIN) < 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance less than 0.0");
        }
        if (varStats.getValue(afwMath::MIN) == 0.0) {
            throw LSST_EXCEPT(pexExcept::Exception, 
                              "Error: variance equals 0.0, cannot inverse variance weight");
        }

        lsst::afw::math::KernelList basisList = 
            boost::shared_dynamic_cast<afwMath::LinearCombinationKernel>(this->_kernel)->getKernelList();

        unsigned int const nKernelParameters     = basisList.size();
        unsigned int const nBackgroundParameters = this->_fitForBackground ? 1 : 0;
        unsigned int const nParameters           = nKernelParameters + nBackgroundParameters;

        std::vector<boost::shared_ptr<afwMath::Kernel> >::const_iterator kiter = basisList.begin();

        /* 
           NOTE : If we are using these views in Afw's Image space, we need to
           make sure and compensate for the XY0 of the image:

           afwGeom::Box2I fullBBox = imageToConvolve.getBBox(afwImage::PARENT);
           int maskStartCol = maskBox.getMinX();
           int maskStartRow = maskBox.getMinY();
           int maskEndCol   = maskBox.getMaxX();
           int maskEndRow   = maskBox.getMaxY();

        
          If we are going to be doing the slicing in Eigen matrices derived from
          the images, ignore the XY0.

           afwGeom::Box2I fullBBox = imageToConvolve.getBBox(afwImage::LOCAL);

           int maskStartCol = maskBox.getMinX() - imageToConvolve.getX0();
           int maskStartRow = maskBox.getMinY() - imageToConvolve.getY0();
           int maskEndCol   = maskBox.getMaxX() - imageToConvolve.getX0();
           int maskEndRow   = maskBox.getMaxY() - imageToConvolve.getY0();

        */

                                                                         
        afwGeom::Box2I shrunkBBox = (*kiter)->shrinkBBox(imageToConvolve.getBBox(afwImage::PARENT));

        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Limits of good pixels after convolution: %d,%d -> %d,%d", 
                          shrunkBBox.getMinX(), shrunkBBox.getMinY(), 
                          shrunkBBox.getMaxX(), shrunkBBox.getMaxY());

        unsigned int const startCol = shrunkBBox.getMinX();
        unsigned int const startRow = shrunkBBox.getMinY();
        unsigned int const endCol   = shrunkBBox.getMaxX();
        unsigned int const endRow   = shrunkBBox.getMaxY();

        /* NOTE: no endCol/endRow += 1 for block slicing, since we are doing the
           slicing using Afw, not Eigen

           Eigen arrays have index 0,0 in the upper right, while LSST images
           have 0,0 in the lower left.  The y-axis is flipped.  When translating
           Images to Eigen matrices in ipDiffim::imageToEigenMatrix this is
           accounted for.  However, we still need to be aware of this fact if
           addressing subregions of an Eigen matrix.  This is why the slicing is
           done in Afw, its cleaner.
           
           Please see examples/maskedKernel.cc for elaboration on some of the
           tests done to make sure this indexing gets done right.

        */


        /* Inner limits; of mask box */
        int maskStartCol = maskBox.getMinX();
        int maskStartRow = maskBox.getMinY();
        int maskEndCol   = maskBox.getMaxX();
        int maskEndRow   = maskBox.getMaxY();

        /* 

        |---------------------------|
        |      Kernel Boundary      |
        |  |---------------------|  |
        |  |         Top         |  |
        |  |......_________......|  |
        |  |      |       |      |  |
        |  |  L   |  Box  |  R   |  |
        |  |      |       |      |  |
        |  |......---------......|  |
        |  |        Bottom       |  |
        |  |---------------------|  |
        |                           |
        |---------------------------|
       
        4 regions we want to extract from the pixels: top bottom left right

        */
        afwGeom::Box2I tBox = afwGeom::Box2I(afwGeom::Point2I(startCol, maskEndRow + 1),
                                             afwGeom::Point2I(endCol, endRow));
        
        afwGeom::Box2I bBox = afwGeom::Box2I(afwGeom::Point2I(startCol, startRow),
                                             afwGeom::Point2I(endCol, maskStartRow - 1));
        
        afwGeom::Box2I lBox = afwGeom::Box2I(afwGeom::Point2I(startCol, maskStartRow),
                                             afwGeom::Point2I(maskStartCol - 1, maskEndRow));
        
        afwGeom::Box2I rBox = afwGeom::Box2I(afwGeom::Point2I(maskEndCol + 1, maskStartRow),
                                             afwGeom::Point2I(endCol, maskEndRow));

        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Upper good pixel region: %d,%d -> %d,%d", 
                          tBox.getMinX(), tBox.getMinY(), tBox.getMaxX(), tBox.getMaxY());
        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Bottom good pixel region: %d,%d -> %d,%d", 
                          bBox.getMinX(), bBox.getMinY(), bBox.getMaxX(), bBox.getMaxY());
        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Left good pixel region: %d,%d -> %d,%d", 
                          lBox.getMinX(), lBox.getMinY(), lBox.getMaxX(), lBox.getMaxY());
        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Right good pixel region: %d,%d -> %d,%d", 
                          rBox.getMinX(), rBox.getMinY(), rBox.getMaxX(), rBox.getMaxY());

        std::vector<afwGeom::Box2I> boxArray;
        boxArray.push_back(tBox);
        boxArray.push_back(bBox);
        boxArray.push_back(lBox);
        boxArray.push_back(rBox);        

        int totalSize = tBox.getWidth() * tBox.getHeight();
        totalSize    += bBox.getWidth() * bBox.getHeight();
        totalSize    += lBox.getWidth() * lBox.getHeight();
        totalSize    += rBox.getWidth() * rBox.getHeight();

        Eigen::MatrixXd eigenToConvolve(totalSize, 1);
        Eigen::MatrixXd eigenToNotConvolve(totalSize, 1);
        Eigen::MatrixXd eigeniVariance(totalSize, 1);
        eigenToConvolve.setZero();
        eigenToNotConvolve.setZero();
        eigeniVariance.setZero();
        
        boost::timer t;
        t.restart();
 
        int nTerms = 0;
        typename std::vector<afwGeom::Box2I>::iterator biter = boxArray.begin();
        for (; biter != boxArray.end(); ++biter) {
            int area = (*biter).getWidth() * (*biter).getHeight();

            afwImage::Image<InputT> siToConvolve(imageToConvolve, *biter, afwImage::PARENT);
            afwImage::Image<InputT> siToNotConvolve(imageToNotConvolve, *biter, afwImage::PARENT);
            afwImage::Image<InputT> sVarEstimate(varianceEstimate, *biter, afwImage::PARENT);

            Eigen::MatrixXd eToConvolve = imageToEigenMatrix(siToConvolve);
            Eigen::MatrixXd eToNotConvolve = imageToEigenMatrix(siToNotConvolve);
            Eigen::MatrixXd eiVarEstimate = imageToEigenMatrix(sVarEstimate).array().inverse().matrix();
            
            eToConvolve.resize(area, 1);
            eToNotConvolve.resize(area, 1);
            eiVarEstimate.resize(area, 1);

            eigenToConvolve.block(nTerms, 0, area, 1) = eToConvolve;
            eigenToNotConvolve.block(nTerms, 0, area, 1) = eToNotConvolve;
            eigeniVariance.block(nTerms, 0, area, 1) = eiVarEstimate;

            nTerms += area;
        }

        afwImage::Image<InputT> cimage(imageToConvolve.getDimensions());

        std::vector<boost::shared_ptr<Eigen::MatrixXd> > convolvedEigenList(nKernelParameters);
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiter = 
            convolvedEigenList.begin();
        /* Create C_i in the formalism of Alard & Lupton */
        for (kiter = basisList.begin(); kiter != basisList.end(); ++kiter, ++eiter) {
            afwMath::convolve(cimage, imageToConvolve, **kiter, false); /* cimage stores convolved image */
            Eigen::MatrixXd cMat(totalSize, 1);
            cMat.setZero();

            int nTerms = 0;
            typename std::vector<afwGeom::Box2I>::iterator biter = boxArray.begin();
            for (; biter != boxArray.end(); ++biter) {
                int area = (*biter).getWidth() * (*biter).getHeight();

                afwImage::Image<InputT> csubimage(cimage, *biter, afwImage::PARENT);
                Eigen::MatrixXd esubimage = imageToEigenMatrix(csubimage);
                esubimage.resize(area, 1);
                cMat.block(nTerms, 0, area, 1) = esubimage;

                nTerms += area;
            }
            
            boost::shared_ptr<Eigen::MatrixXd> cMatPtr (new Eigen::MatrixXd(cMat));
            *eiter = cMatPtr;
            
        } 
        
        double time = t.elapsed();
        pexLog::TTrace<5>("lsst.ip.diffim.MaskedKernelSolution.build", 
                          "Total compute time to do basis convolutions : %.2f s", time);
        t.restart();

        /* 
           Load matrix with all values from convolvedEigenList : all images
           (eigeniVariance, convolvedEigenList) must be the same size
        */
        Eigen::MatrixXd cMat(eigenToConvolve.col(0).size(), nParameters);
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterj = 
            convolvedEigenList.begin();
        typename std::vector<boost::shared_ptr<Eigen::MatrixXd> >::iterator eiterE = 
            convolvedEigenList.end();
        for (unsigned int kidxj = 0; eiterj != eiterE; eiterj++, kidxj++) {
            cMat.col(kidxj) = (*eiterj)->col(0);
        }
        /* Treat the last "image" as all 1's to do the background calculation. */
        if (this->_fitForBackground)
            cMat.col(nParameters-1).fill(1.);

        this->_cMat.reset(new Eigen::MatrixXd(cMat));
        this->_ivVec.reset(new Eigen::VectorXd(eigeniVariance.col(0)));
        this->_iVec.reset(new Eigen::VectorXd(eigenToNotConvolve.col(0)));

        /* Make these outside of solve() so I can check condition number */
        this->_mMat.reset(
            new Eigen::MatrixXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_cMat))
            );
        this->_bVec.reset(
            new Eigen::VectorXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_iVec))
            );
        
    }
    /*******************************************************************************************************/


    template <typename InputT>
    RegularizedKernelSolution<InputT>::RegularizedKernelSolution(
        lsst::afw::math::KernelList const& basisList,
        bool fitForBackground,
        boost::shared_ptr<Eigen::MatrixXd> hMat,
        lsst::pex::policy::Policy policy
        ) 
        :
        StaticKernelSolution<InputT>(basisList, fitForBackground),
        _hMat(hMat),
        _policy(policy)
    {};

    template <typename InputT>
    double RegularizedKernelSolution<InputT>::estimateRisk(double maxCond) {
        Eigen::MatrixXd vMat      = (this->_cMat)->jacobiSvd().matrixV();
        Eigen::MatrixXd vMatvMatT = vMat * vMat.transpose();

        /* Find pseudo inverse of mMat, which may be ill conditioned */
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVecValues(*(this->_mMat));
        Eigen::MatrixXd const& rMat = eVecValues.eigenvectors();
        Eigen::VectorXd eValues = eVecValues.eigenvalues();
        double eMax = eValues.maxCoeff();
        for (int i = 0; i != eValues.rows(); ++i) {
            if (eValues(i) == 0.0) {
                eValues(i) = 0.0;
            }
            else if ((eMax / eValues(i)) > maxCond) {
                pexLog::TTrace<5>("lsst.ip.diffim.RegularizedKernelSolution.estimateRisk", 
                                  "Truncating eValue %d; %.5e / %.5e = %.5e vs. %.5e",
                                  i, eMax, eValues(i), eMax / eValues(i), maxCond);
                eValues(i) = 0.0;
            }
            else {
                eValues(i) = 1.0 / eValues(i);
            }
        }
        Eigen::MatrixXd mInv    = rMat * eValues.asDiagonal() * rMat.transpose();

        std::vector<double> lambdas = _createLambdaSteps();
        std::vector<double> risks;
        for (unsigned int i = 0; i < lambdas.size(); i++) {
            double l = lambdas[i];
            Eigen::MatrixXd mLambda = *(this->_mMat) + l * (*_hMat);
            
            try {
                KernelSolution::solve(mLambda, *(this->_bVec));
            } catch (pexExcept::Exception &e) {
                LSST_EXCEPT_ADD(e, "Unable to solve regularized kernel matrix");
                throw e;
            }
            Eigen::VectorXd term1 = (this->_aVec->transpose() * vMatvMatT * *(this->_aVec));
            if (term1.size() != 1)
                throw LSST_EXCEPT(pexExcept::Exception, "Matrix size mismatch");

            double term2a = (vMatvMatT * mLambda.inverse()).trace();

            Eigen::VectorXd term2b = (this->_aVec->transpose() * (mInv * *(this->_bVec)));
            if (term2b.size() != 1)
                throw LSST_EXCEPT(pexExcept::Exception, "Matrix size mismatch");

            double risk   = term1(0) + 2 * (term2a - term2b(0));
            pexLog::TTrace<6>("lsst.ip.diffim.RegularizedKernelSolution.estimateRisk", 
                              "Lambda = %.3f, Risk = %.5e", 
                              l, risk);
            pexLog::TTrace<7>("lsst.ip.diffim.RegularizedKernelSolution.estimateRisk", 
                              "%.5e + 2 * (%.5e - %.5e)", 
                              term1(0), term2a, term2b(0));
            risks.push_back(risk);
        }
        std::vector<double>::iterator it = min_element(risks.begin(), risks.end());
        int index = distance(risks.begin(), it);
        pexLog::TTrace<5>("lsst.ip.diffim.RegularizedKernelSolution.estimateRisk", 
                          "Minimum Risk = %.3e at lambda = %.3e", risks[index], lambdas[index]);

        return lambdas[index];
        
    }

    template <typename InputT>
    boost::shared_ptr<Eigen::MatrixXd> RegularizedKernelSolution<InputT>::getM(bool includeHmat) {
        if (includeHmat == true) {
            return (boost::shared_ptr<Eigen::MatrixXd>(
                        new Eigen::MatrixXd(*(this->_mMat) + _lambda * (*_hMat))
                        ));
        }
        else {
            return this->_mMat;
        }
    }

    template <typename InputT>
    void RegularizedKernelSolution<InputT>::solve() {

        pexLog::TTrace<5>("lsst.ip.diffim.RegularizedKernelSolution.solve", 
                          "cMat is %d x %d; vVec is %d; iVec is %d; hMat is %d x %d", 
                          (*this->_cMat).rows(), (*this->_cMat).cols(), (*this->_ivVec).size(), 
                          (*this->_iVec).size(), (*_hMat).rows(), (*_hMat).cols());

        if (DEBUG_MATRIX2) {
            std::cout << "ID: " << (this->_id) << std::endl;
            std::cout << "C:" << std::endl;
            std::cout << (*this->_cMat) << std::endl;
            std::cout << "Sigma^{-1}:" << std::endl;
            std::cout << Eigen::MatrixXd((*this->_ivVec).asDiagonal()) << std::endl;
            std::cout << "Y:" << std::endl;
            std::cout << (*this->_iVec) << std::endl;
            std::cout << "H:" << std::endl;
            std::cout << (*_hMat) << std::endl;
        }


        this->_mMat.reset(
            new Eigen::MatrixXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_cMat))
            );
        this->_bVec.reset(
            new Eigen::VectorXd(this->_cMat->transpose() * this->_ivVec->asDiagonal() * *(this->_iVec))
            );
        
        std::string lambdaType = _policy.getString("lambdaType");        
        double lambdaValue     = _policy.getDouble("lambdaValue");
        
        /* See N.R. 18.5
           
           Matrix equation to solve is Y = C a + N
           Y   = vectorized version of I (I = image to not convolve)
           C_i = K_i (x) R (R = image to convolve)
           a   = kernel coefficients
           N   = noise
           
           If we reweight everything by the inverse square root of the noise
           covariance, we get a linear model with the identity matrix for
           the noise.  The problem can then be solved using least squares,
           with the normal equations
           
              C^T Y = C^T C a

           or 
           
              b = M a

           with 

              b = C^T Y
              M = C^T C
              a = (C^T C)^{-1} C^T Y


           We can regularize the least square problem
  
              Y = C a + lambda a^T H a       (+ N, which can be ignored)

           or the normal equations

              (C^T C + lambda H) a = C^T Y
 
              
           The solution to the regularization of the least squares problem is

              a = (C^T C + lambda H)^{-1} C^T Y

           The approximation to Y is 

              C a = C (C^T C + lambda H)^{-1} C^T Y
 
           with smoothing matrix 

              S = C (C^T C + lambda H)^{-1} C^T 

        */
        
        if (lambdaType == "absolute") {
            _lambda = lambdaValue;
        }
        else if (lambdaType ==  "relative") {
            _lambda  = this->_mMat->trace() / this->_hMat->trace();
            _lambda *= lambdaValue;
        }
        else if (lambdaType ==  "minimizeBiasedRisk") {
            double tol = _policy.getDouble("maxConditionNumber");
            _lambda = estimateRisk(tol);
        }
        else if (lambdaType ==  "minimizeUnbiasedRisk") {
            _lambda = estimateRisk(std::numeric_limits<double>::max());
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "lambdaType in Policy not recognized");
        }
        
        pexLog::TTrace<5>("lsst.ip.diffim.KernelCandidate.build", 
                          "Applying kernel regularization with lambda = %.2e", _lambda);
        
        
        try {
            KernelSolution::solve(*(this->_mMat) + _lambda * (*_hMat), *(this->_bVec));
        } catch (pexExcept::Exception &e) {
            LSST_EXCEPT_ADD(e, "Unable to solve static kernel matrix");
            throw e;
        }
        /* Turn matrices into _kernel and _background */
        StaticKernelSolution<InputT>::_setKernel();
    }

    template <typename InputT>
    std::vector<double> RegularizedKernelSolution<InputT>::_createLambdaSteps() {
        std::vector<double> lambdas;

        std::string lambdaStepType = _policy.getString("lambdaStepType");    
        if (lambdaStepType == "linear") {
            double lambdaLinMin   = _policy.getDouble("lambdaLinMin");        
            double lambdaLinMax   = _policy.getDouble("lambdaLinMax");
            double lambdaLinStep  = _policy.getDouble("lambdaLinStep");
            for (double l = lambdaLinMin; l <= lambdaLinMax; l += lambdaLinStep) {
                lambdas.push_back(l);
            }
        }
        else if (lambdaStepType == "log") {
            double lambdaLogMin   = _policy.getDouble("lambdaLogMin");        
            double lambdaLogMax   = _policy.getDouble("lambdaLogMax");
            double lambdaLogStep  = _policy.getDouble("lambdaLogStep");
            for (double l = lambdaLogMin; l <= lambdaLogMax; l += lambdaLogStep) {
                lambdas.push_back(pow(10, l));
            }
        }
        else {
            throw LSST_EXCEPT(pexExcept::Exception, "lambdaStepType in Policy not recognized");
        }
        return lambdas;
    }

    /*******************************************************************************************************/

    SpatialKernelSolution::SpatialKernelSolution(
        lsst::afw::math::KernelList const& basisList,
        lsst::afw::math::Kernel::SpatialFunctionPtr spatialKernelFunction,
        lsst::afw::math::Kernel::SpatialFunctionPtr background,
        lsst::pex::policy::Policy policy
        ) :
        KernelSolution(),
        _spatialKernelFunction(spatialKernelFunction),
        _constantFirstTerm(false),
        _kernel(),
        _background(background),
        _kSum(0.0),
        _policy(policy),
        _nbases(0),
        _nkt(0),
        _nbt(0),
        _nt(0) {

        bool isAlardLupton    = _policy.getString("kernelBasisSet") == "alard-lupton";
        bool usePca           = _policy.getBool("usePcaForSpatialKernel");
        if (isAlardLupton || usePca) {
            _constantFirstTerm = true;
        }
        this->_fitForBackground = _policy.getBool("fitForBackground");

        _nbases = basisList.size();
        _nkt = _spatialKernelFunction->getParameters().size();
        _nbt = _fitForBackground ? _background->getParameters().size() : 0;
        _nt  = 0;
        if (_constantFirstTerm) {
            _nt = (_nbases - 1) * _nkt + 1 + _nbt;
        } else {
            _nt = _nbases * _nkt + _nbt;
        }
        
        boost::shared_ptr<Eigen::MatrixXd> mMat (new Eigen::MatrixXd(_nt, _nt));
        boost::shared_ptr<Eigen::VectorXd> bVec (new Eigen::VectorXd(_nt));
        (*mMat).setZero();
        (*bVec).setZero();

        this->_mMat = mMat;
        this->_bVec = bVec;

        _kernel = afwMath::LinearCombinationKernel::Ptr(
            new afwMath::LinearCombinationKernel(basisList, *_spatialKernelFunction)
            );

        pexLog::TTrace<5>("lsst.ip.diffim.SpatialKernelSolution", 
                          "Initializing with size %d %d %d and constant first term = %s",
                          _nkt, _nbt, _nt,
                          _constantFirstTerm ? "true" : "false");
        
    }

    void SpatialKernelSolution::addConstraint(float xCenter, float yCenter,
                                              boost::shared_ptr<Eigen::MatrixXd> qMat,
                                              boost::shared_ptr<Eigen::VectorXd> wVec) {
        
        pexLog::TTrace<8>("lsst.ip.diffim.SpatialKernelSolution.addConstraint", 
                          "Adding candidate at %f, %f", xCenter, yCenter);

        /* Calculate P matrices */
        /* Pure kernel terms */
        Eigen::VectorXd pK(_nkt);
        std::vector<double> paramsK = _spatialKernelFunction->getParameters();
        for (int idx = 0; idx < _nkt; idx++) { paramsK[idx] = 0.0; }
        for (int idx = 0; idx < _nkt; idx++) {
            paramsK[idx] = 1.0;
            _spatialKernelFunction->setParameters(paramsK);
            pK(idx) = (*_spatialKernelFunction)(xCenter, yCenter); /* Assume things don't vary over stamp */
            paramsK[idx] = 0.0;
        }
        Eigen::MatrixXd pKpKt = (pK * pK.transpose());
        
        Eigen::VectorXd pB;
        Eigen::MatrixXd pBpBt;
        Eigen::MatrixXd pKpBt;
        if (_fitForBackground) {
            pB = Eigen::VectorXd(_nbt);

            /* Pure background terms */
            std::vector<double> paramsB = _background->getParameters();
            for (int idx = 0; idx < _nbt; idx++) { paramsB[idx] = 0.0; }
            for (int idx = 0; idx < _nbt; idx++) {
                paramsB[idx] = 1.0;
                _background->setParameters(paramsB);
                pB(idx) = (*_background)(xCenter, yCenter);       /* Assume things don't vary over stamp */
                paramsB[idx] = 0.0;
            }
            pBpBt = (pB * pB.transpose());
            
            /* Cross terms */
            pKpBt = (pK * pB.transpose());
        }
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial weights" << std::endl;
            std::cout << "pKpKt " << pKpKt << std::endl;
            if (_fitForBackground) {
                std::cout << "pBpBt " << pBpBt << std::endl;
                std::cout << "pKpBt " << pKpBt << std::endl;
            }
        }
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial matrix inputs" << std::endl;
            std::cout << "M " << (*qMat) << std::endl;
            std::cout << "B " << (*wVec) << std::endl;
        }

        /* first index to start the spatial blocks; default=0 for non-constant first term */
        int m0 = 0; 
        /* how many rows/cols to adjust the matrices/vectors; default=0 for non-constant first term */
        int dm = 0; 
        /* where to start the background terms; this is always true */
        int mb = _nt - _nbt;
        
        if (_constantFirstTerm) {
            m0 = 1;       /* we need to manually fill in the first (non-spatial) terms below */
            dm = _nkt-1;  /* need to shift terms due to lack of spatial variation in first term */
            
            (*_mMat)(0, 0) += (*qMat)(0,0);
            for(int m2 = 1; m2 < _nbases; m2++)  {
                (*_mMat).block(0, m2*_nkt-dm, 1, _nkt) += (*qMat)(0,m2) * pK.transpose();
            }
            (*_bVec)(0) += (*wVec)(0);
            
            if (_fitForBackground) {
                (*_mMat).block(0, mb, 1, _nbt) += (*qMat)(0,_nbases) * pB.transpose();
            }
        }
        
        /* Fill in the spatial blocks */
        for(int m1 = m0; m1 < _nbases; m1++)  {
            /* Diagonal kernel-kernel term; only use upper triangular part of pKpKt */
            (*_mMat).block(m1*_nkt-dm, m1*_nkt-dm, _nkt, _nkt) += 
                (pKpKt * (*qMat)(m1,m1)).triangularView<Eigen::Upper>();
            
            /* Kernel-kernel terms */
            for(int m2 = m1+1; m2 < _nbases; m2++)  {
                (*_mMat).block(m1*_nkt-dm, m2*_nkt-dm, _nkt, _nkt) += (*qMat)(m1,m2) * pKpKt;
            }

            if (_fitForBackground) {
                /* Kernel cross terms with background */
                (*_mMat).block(m1*_nkt-dm, mb, _nkt, _nbt) += (*qMat)(m1,_nbases) * pKpBt;
            }
            
            /* B vector */
            (*_bVec).segment(m1*_nkt-dm, _nkt) += (*wVec)(m1) * pK;
        }
        
        if (_fitForBackground) {
            /* Background-background terms only */
            (*_mMat).block(mb, mb, _nbt, _nbt) +=  
                (pBpBt * (*qMat)(_nbases,_nbases)).triangularView<Eigen::Upper>();
            (*_bVec).segment(mb, _nbt)         += (*wVec)(_nbases) * pB;
        }
        
        if (DEBUG_MATRIX) {
            std::cout << "Spatial matrix outputs" << std::endl;
            std::cout << "mMat " << (*_mMat) << std::endl;
            std::cout << "bVec " << (*_bVec) << std::endl;
        }

    }

    KernelSolution::ImageT::Ptr SpatialKernelSolution::makeKernelImage() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return image");
        }
        ImageT::Ptr image(new ImageT(_kernel->getDimensions()));
        (void)_kernel->computeImage(*image, false);              
        return image;
    }

    void SpatialKernelSolution::solve() {
        /* Fill in the other half of mMat */
        for (int i = 0; i < _nt; i++) {
            for (int j = i+1; j < _nt; j++) {
                (*_mMat)(j,i) = (*_mMat)(i,j);
            }
        }

        try {
            KernelSolution::solve();
        } catch (pexExcept::Exception &e) {
            LSST_EXCEPT_ADD(e, "Unable to solve spatial kernel matrix");
            throw e;
        }
        /* Turn matrices into _kernel and _background */
        _setKernel();
    }

    std::pair<afwMath::LinearCombinationKernel::Ptr, 
            afwMath::Kernel::SpatialFunctionPtr> SpatialKernelSolution::getSolutionPair() {
        if (_solvedBy == KernelSolution::NONE) {
            throw LSST_EXCEPT(pexExcept::Exception, "Kernel not solved; cannot return solution");
        }
        
        return std::make_pair(_kernel, _background);
    }
    
    void SpatialKernelSolution::_setKernel() {
        /* Report on condition number */
        double cNumber = this->getConditionNumber(EIGENVALUE);
        
        if (_nkt == 1) {
            /* Not spatially varying; this fork is a specialization for convolution speed--up */
            
            /* Set the basis coefficients */
            std::vector<double> kCoeffs(_nbases);
            for (int i = 0; i < _nbases; i++) {
                if (std::isnan((*_aVec)(i))) {
                    throw LSST_EXCEPT(
                        pexExcept::Exception, 
                        str(boost::format(
                                "I. Unable to determine spatial kernel solution %d (nan).  Condition number = %.3e") % i % cNumber));
                }
                kCoeffs[i] = (*_aVec)(i);
            }
            lsst::afw::math::KernelList basisList = 
                boost::shared_dynamic_cast<afwMath::LinearCombinationKernel>(_kernel)->getKernelList();
            _kernel.reset(
                new afwMath::LinearCombinationKernel(basisList, kCoeffs)
                );
        }
        else {
            
            /* Set the kernel coefficients */
            std::vector<std::vector<double> > kCoeffs;
            kCoeffs.reserve(_nbases);
            for (int i = 0, idx = 0; i < _nbases; i++) {
                kCoeffs.push_back(std::vector<double>(_nkt));

                /* Deal with the possibility the first term doesn't vary spatially */
                if ((i == 0) && (_constantFirstTerm)) {
                    if (std::isnan((*_aVec)(idx))) {
                        throw LSST_EXCEPT(
                            pexExcept::Exception, 
                            str(boost::format(
                                    "II. Unable to determine spatial kernel solution %d (nan).  Condition number = %.3e") % idx % cNumber));
                    }
                    kCoeffs[i][0] = (*_aVec)(idx++);
                }
                else {
                    for (int j = 0; j < _nkt; j++) {
                        if (std::isnan((*_aVec)(idx))) {
                            throw LSST_EXCEPT(
                                pexExcept::Exception, 
                                str(boost::format(
                                        "III. Unable to determine spatial kernel solution %d (nan).  Condition number = %.3e") % idx % cNumber));
                        }
                        kCoeffs[i][j] = (*_aVec)(idx++);
                    }
                }
            }
            _kernel->setSpatialParameters(kCoeffs);
        }

        /* Set kernel Sum */
        ImageT::Ptr image (new ImageT(_kernel->getDimensions()));
        _kSum  = _kernel->computeImage(*image, false);              

        /* Set the background coefficients */
        std::vector<double> bgCoeffs(_fitForBackground ? _nbt : 1);
        if (_fitForBackground) {
            for (int i = 0; i < _nbt; i++) {
                int idx = _nt - _nbt + i;
                if (std::isnan((*_aVec)(idx))) {
                    throw LSST_EXCEPT(pexExcept::Exception, 
                                      str(boost::format(
                           "Unable to determine spatial background solution %d (nan)") % (idx)));
                }
                bgCoeffs[i] = (*_aVec)(idx);
            }
        }
        else {
            bgCoeffs[0] = 0.;
        }
        _background->setParameters(bgCoeffs);
    }

/***********************************************************************************************************/
//
// Explicit instantiations
//
    typedef float InputT;

    template class StaticKernelSolution<InputT>;
    template class MaskedKernelSolution<InputT>;
    template class RegularizedKernelSolution<InputT>;

}}} // end of namespace lsst::ip::diffim
