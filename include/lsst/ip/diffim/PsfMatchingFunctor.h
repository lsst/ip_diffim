// -*- lsst-c++ -*-
/**
 * @file PsfMatchingFunctor.h
 *
 * @brief Class owning and implementing the core functionality of building a single Psf matching kernel
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_PSFMATCHINGFUNCTOR_H
#define LSST_IP_DIFFIM_PSFMATCHINGFUNCTOR_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/image/MaskedImage.h>

namespace pexPolicy  = lsst::pex::policy; 
namespace afwMath    = lsst::afw::math;
namespace afwImage   = lsst::afw::image;

namespace lsst { namespace ip { namespace diffim {
   
   /**
     * @brief Functor to create PSF Matching Kernel
     *
     * @note This class owns the functionality to make a single difference
     * imaging kernel around one object realized in 2 different images.  If
     * constructed with a regularization matrix, will use it by default.  This
     * creates the M and B vectors that are used to solve for the kernel
     * parameters 'x' as in Mx = B.  This creates a single kernel around a
     * single object, and operates in tandem with the KernelCandidate +
     * BuildSingleKernelVisitor classes for the spatial modeling.
     * 
     * @ingroup ip_diffim
     */
    template <typename PixelT, typename VarT=afwImage::VariancePixel>
    class PsfMatchingFunctor {
    public:
        typedef boost::shared_ptr<PsfMatchingFunctor> Ptr;
        typedef typename afwImage::MaskedImage<PixelT>::xy_locator xy_locator;
        typedef typename afwImage::Image<VarT>::xy_locator         xyi_locator;

        PsfMatchingFunctor(
            afwMath::KernelList const& basisList
            );
        PsfMatchingFunctor(
            afwMath::KernelList const& basisList,
            boost::shared_ptr<Eigen::MatrixXd> const& _hMat
            );
        virtual ~PsfMatchingFunctor() {};

        /* Shallow copy only; shared matrix product uninitialized */
        PsfMatchingFunctor(const PsfMatchingFunctor<PixelT,VarT> &rhs);

        std::pair<boost::shared_ptr<afwMath::Kernel>, double> getSolution();
        std::pair<boost::shared_ptr<afwMath::Kernel>, double> getSolutionUncertainty();
        
        /** Access to least squares info
         */
        std::pair<boost::shared_ptr<Eigen::MatrixXd>, boost::shared_ptr<Eigen::VectorXd> > getAndClearMB();

        /** Access to basis list
         */
        afwMath::KernelList getBasisList() const { return _basisList; }

        /* Create PSF matching kernel */
        void apply(afwImage::Image<PixelT> const& imageToConvolve,
                   afwImage::Image<PixelT> const& imageToNotConvolve,
                   afwImage::Image<VarT>   const& varianceEstimate,
                   pexPolicy::Policy       const& policy
            );

    protected:
        afwMath::KernelList const _basisList;                    ///< List of Kernel basis functions
        boost::shared_ptr<Eigen::MatrixXd> _mMat;                ///< Least squares matrix
        boost::shared_ptr<Eigen::VectorXd> _bVec;                ///< Least squares vector
        boost::shared_ptr<Eigen::VectorXd> _sVec;                ///< Least square solution
        boost::shared_ptr<Eigen::MatrixXd> const _hMat;          ///< Regularization matrix
        bool _initialized;                                       ///< Has been solved for
        bool _regularize;                                        ///< Has a _H matrix
    };
    
    /**
     * @brief Helper method to return a pointer to a PsfMatchingFunctor()
     *
     * @param basisList  Input set of basis kernels to use for Psf matching
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    typename PsfMatchingFunctor<PixelT>::Ptr
    makePsfMatchingFunctor(afwMath::KernelList const& basisList) {
        return typename PsfMatchingFunctor<PixelT>::Ptr(new PsfMatchingFunctor<PixelT>(basisList));
    }

    /**
     * @brief Helper method to return a pointer to a PsfMatchingFunctor() with regularization
     *
     * @param basisList  Input set of basis kernels to use for Psf matching
     * @param H  Regularization matrix (for delta-function bases)
     *
     * @ingroup ip_diffim
     */
    template <typename PixelT>
    typename PsfMatchingFunctor<PixelT>::Ptr
    makePsfMatchingFunctor(afwMath::KernelList const& basisList,
                           boost::shared_ptr<Eigen::MatrixXd> const _hMat) {
        return typename PsfMatchingFunctor<PixelT>::Ptr(new PsfMatchingFunctor<PixelT>(basisList, _hMat));
    }

}}}

#endif
