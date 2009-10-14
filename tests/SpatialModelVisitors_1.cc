#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SpatialModelVisitors

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/pex/policy/Policy.h"
#include "lsst/ip/diffim/SpatialModelVisitors.h"

using lsst::pex::policy::Policy;

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
using namespace lsst::ip::diffim;

typedef float PixelT;
typedef afwImage::Image<afwMath::Kernel::Pixel> ImageT;

BOOST_AUTO_TEST_CASE(KernelSumVisitor) {
    Policy::Ptr policy(new Policy);
    policy->set("kernelSumClipping", false);
    policy->set("maxKsumSigma", 3.0);

    image::MaskedImage<PixelT>::Ptr mimg1(
        new image::MaskedImage<PixelT>(100,100)
        );
    image::MaskedImage<PixelT>::Ptr mimg2(
        new image::MaskedImage<PixelT>(100,100)
        );

    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand1(new KernelCandidate<PixelT>(0., 0., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand2(new KernelCandidate<PixelT>(0., 0., mimg1, mimg2));
    
    ImageT img1(10,10);
    img1 = 0;
    *img1.at(4, 4) = 1;

    ImageT img2(10,10);
    img2 = 0;
    *img2.at(4, 4) = 100;
    
    afwMath::Kernel::Ptr k1(new afwMath::FixedKernel(img1));
    afwMath::Kernel::Ptr k2(new afwMath::FixedKernel(img2));

    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->setKernel(k1);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand2))->setKernel(k2);
    
    /* Basic sum */
    {
        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processKsumDistribution();
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumMean(),  101*0.5);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumNpts(),  2);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getNRejected(), 0);   
    }
    
    /* There is sigma clipping in the stats, but kernelSumClipping = false */
    {
        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processKsumDistribution();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumMean(),  1.0);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumNpts(),  4);   /* cand2 sigma clipped out in the stats */
        BOOST_CHECK_EQUAL(kernelSumVisitor.getNRejected(), 0);   /* but not from the cells */
        BOOST_CHECK_EQUAL(cand2->getStatus(), afwMath::SpatialCellCandidate::UNKNOWN);
    }

    /* There is sigma clipping in the stats, but kernelSumClipping = false */
    {
        policy->set("kernelSumClipping", true);
        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processKsumDistribution();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumMean(),  1.0);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumNpts(),  4);   /* cand2 sigma clipped out in the stats */
        BOOST_CHECK_EQUAL(kernelSumVisitor.getNRejected(), 1);   /* and from the cells */
        BOOST_CHECK_EQUAL(cand2->getStatus(), afwMath::SpatialCellCandidate::BAD);
    }
}

#if 0
BOOST_AUTO_TEST_CASE(SetPcaImageVisitor) {
    afwImage::ImagePca<ImageT> imagePca;
    detail::SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
}

BOOST_AUTO_TEST_CASE(BuildSingleKernelVisitor) {
    detail::BuildSingleKernelVisitor<PixelT> singleKernelFitterPca(kFunctorPca, policy);
}

BOOST_AUTO_TEST_CASE(AssessSpatialKernelVisitor) {
    detail::AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, spatialBackground, policy);
}

BOOST_AUTO_TEST_CASE(BuildSpatialKernelVisitor) {
    detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse, spatialKernelOrder, spatialBgOrder, policy);
}
#endif
