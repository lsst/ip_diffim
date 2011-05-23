// -*- lsst-c++ -*-
/**
 * @file KernelSumVisitor.h
 *
 * @brief Declaration of KernelSumVisitor 
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_KERNELSUMVISITOR_H
#define LSST_IP_DIFFIM_KERNELSUMVISITOR_H

#include "boost/shared_ptr.hpp"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/pex/policy/Policy.h"

namespace lsst { 
namespace ip { 
namespace diffim { 
namespace detail {

    template<typename PixelT>
    class KernelSumVisitor : public lsst::afw::math::CandidateVisitor {
    public:
        typedef boost::shared_ptr<KernelSumVisitor<PixelT> > Ptr;

        enum Mode {AGGREGATE = 0, REJECT = 1};
        
        KernelSumVisitor(lsst::pex::policy::Policy const& policy);
        virtual ~KernelSumVisitor() {};
        
        void setMode(Mode mode) {_mode = mode;}
        int    getNRejected() {return _nRejected;}
        double getkSumMean()  {return _kSumMean;}
        double getkSumStd()   {return _kSumStd;}
        double getdkSumMax()  {return _dkSumMax;}
        int    getkSumNpts()  {return _kSumNpts;}
        
        void resetKernelSum();
        void processCandidate(lsst::afw::math::SpatialCellCandidate *candidate);
        void processKsumDistribution();
        
    private:
        Mode _mode;                  ///< Processing mode; AGGREGATE or REJECT
        std::vector<double> _kSums;  ///< List of all candidate kernel sums
        double _kSumMean;            ///< Clipped mean of the kernel sums
        double _kSumStd;             ///< Clipped standard deviation of kernel sums
        double _dkSumMax;            ///< Maximum acceptable deviation from mean sum
        int    _kSumNpts;            ///< Number of points used in the statistics
        int    _nRejected;           ///< Number of candidates rejected during processCandidate()
        lsst::pex::policy::Policy _policy;   ///< Policy controlling behavior
    };    
    
    template<typename PixelT>
    boost::shared_ptr<KernelSumVisitor<PixelT> >
    makeKernelSumVisitor(lsst::pex::policy::Policy const& policy) {
        return typename KernelSumVisitor<PixelT>::Ptr(new KernelSumVisitor<PixelT>(policy));
    }

}}}} // end of namespace lsst::ip::diffim::detail

#endif
