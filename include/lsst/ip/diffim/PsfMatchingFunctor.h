// -*- lsst-c++ -*-
/**
 * @file PsfMatchingFunctor.h
 *
 * @brief Defines the class PsfMatchingFunctor
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup ip_diffim
 */

#ifndef LSST_IP_DIFFIM_PSFMATCHINGFUNCTOR_H
#define LSST_IP_DIFFIM_PSFMATCHINGFUNCTOR_H

#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <boost/timer.hpp> 

#include <boost/shared_ptr.hpp>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>

#include <lsst/pex/policy/Policy.h>
#include <lsst/afw/math.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/math/Function.h>
#include <lsst/afw/math/KernelFunctions.h>
#include <lsst/afw/math/ConvolveImage.h>
#include <lsst/afw/image/Mask.h>
#include <lsst/afw/image/Image.h>
#include <lsst/afw/detection/Footprint.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/pex/logging/Log.h>


namespace lsst{
namespace ip {
namespace diffim{

/** Functor to create PSF Matching Kernel
 *
 * @ingroup diffim
 * 
 */
template <typename ImageT, typename VarT=lsst::afw::image::VariancePixel>
class PsfMatchingFunctor {
public:
    typedef boost::shared_ptr<PsfMatchingFunctor> Ptr;
    typedef lsst::afw::math::Kernel Kernel
    typedef lsst::afw::math::KernelList<Kernel> BasisList;
    typedef lsst::afw::image::Image<ImageT> Image;
    typedef lsst::afw::image::Image<VarT> Variance;
    typedef std::vector<Image::Ptr> ImageList;
    typedef std::vector<Kernel::PtrT> kernelList;
    typedef typename lsst::afw::image::Image<ImageT>::iterator ImgIterator;
    typedef typename lsst::afw::image::Image<VarT>::iterator  VarIterator;

    PsfMatchingFunctor(BasisList const& basisList);
    virtual ~PsfMatchingFunctor() {};

    /** Return background value
     */
    double getBackground() const { return _background; }

    /** Return uncertainty on background value
     */
    double getBackgroundError() const { return _backgroundError; }

    /** Return PSF matching kernel
     */
    KernelPtr getKernel()  const { return _kernel; }

    /** Return uncertainty on matching kernel, as kernel itself
     */
    KernelPtr getKernelError() const { return _kernelError; }

    /** Reset protected class members
     */
    void reset();

    /* Create PSF matching kernel */
    void apply(Image const& imageToConvolve,
               Image const& imageToNotConvolve,
               Variance const& varianceEstimate,
               lsst::pex::policy::Policy const& policy);

protected:
    BasisList _basisList;   ///< List of Kernel basis functions
    double _background; ///< Differenaitl background estimate
    double _backgroundError;///< Uncertainty on background
    KernelPtr _kernel;  ///< PSF matching kernel
    KernelPtr _kernelError; ///< Uncertainty on kernel
};

}}} //end namespace lsst::ip::diffim

#endif
