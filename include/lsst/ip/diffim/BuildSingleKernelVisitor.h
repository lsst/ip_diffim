// -*- lsst-c++ -*-
/**
 * @file BuildSingleKernelVisitor.h
 *
 * @brief Declaration of BuildSingleKernelVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_BUILDSINGLEKERNELVISITOR_H
#define LSST_IP_DIFFIM_BUILDSINGLEKERNELVISITOR_H

#include "boost/shared_ptr.hpp"

#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

#include "lsst/pex/policy/Policy.h"

#include "lsst/ip/diffim/ImageStatistics.h"

namespace lsst { 
namespace ip { 
namespace diffim { 
namespace detail {

    template<typename PixelT>
    class BuildSingleKernelVisitor : public lsst::afw::math::CandidateVisitor {
        typedef lsst::afw::image::MaskedImage<PixelT> MaskedImageT;
    public:
        typedef boost::shared_ptr<BuildSingleKernelVisitor<PixelT> > Ptr;

        BuildSingleKernelVisitor(
            lsst::afw::math::KernelList const& basisList,
            lsst::pex::policy::Policy const& policy  
            );
        BuildSingleKernelVisitor(
            lsst::afw::math::KernelList const& basisList,
            lsst::pex::policy::Policy const& policy, 
            boost::shared_ptr<Eigen::MatrixXd> hMat  
            );
        virtual ~BuildSingleKernelVisitor() {};
        
        /* 
           Don't reprocess candidate if its already been build.  The use
           case for this functionality is : when iterating over all Cells
           and rejecting bad Kernels, we need to re-visit *all* Cells to
           build the next candidate in the list.  Without this flag we would
           unncessarily re-build all the good Kernels.
        */
        void setSkipBuilt(bool skip)      {_skipBuilt = skip;}
        
        int getNRejected()    {return _nRejected;}
        int getNProcessed()   {return _nProcessed;}
        void reset()          {_nRejected = 0; _nProcessed = 0;}
        
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);

    private:
        lsst::afw::math::KernelList const _basisList; ///< Basis set
        lsst::pex::policy::Policy _policy;            ///< Policy controlling behavior
        boost::shared_ptr<Eigen::MatrixXd> _hMat;     ///< Regularization matrix
        ImageStatistics<PixelT> _imstats;     ///< To calculate statistics of difference image
        bool _skipBuilt;                      ///< Skip over built candidates during processCandidate()
        int _nRejected;                       ///< Number of candidates rejected during processCandidate()
        int _nProcessed;                      ///< Number of candidates processed during processCandidate()
        bool _useRegularization;              ///< Regularize if delta function basis

        bool _useCoreStats;                   ///< Extracted from _policy
        int _coreRadius;                      ///< Extracted from _policy
    };
    
    template<typename PixelT>
    boost::shared_ptr<BuildSingleKernelVisitor<PixelT> >
    makeBuildSingleKernelVisitor(
        lsst::afw::math::KernelList const& basisList,
        lsst::pex::policy::Policy const& policy
        ) {

        return typename BuildSingleKernelVisitor<PixelT>::Ptr(
            new BuildSingleKernelVisitor<PixelT>(basisList, policy)
            );
    }

    template<typename PixelT>
    boost::shared_ptr<BuildSingleKernelVisitor<PixelT> >
    makeBuildSingleKernelVisitor(
        lsst::afw::math::KernelList const& basisList,
        lsst::pex::policy::Policy const& policy,
        boost::shared_ptr<Eigen::MatrixXd> hMat  
        ) {

        return typename BuildSingleKernelVisitor<PixelT>::Ptr(
            new BuildSingleKernelVisitor<PixelT>(basisList, policy, hMat)
            );
    }

}}}} // end of namespace lsst::ip::diffim::detail

#endif
