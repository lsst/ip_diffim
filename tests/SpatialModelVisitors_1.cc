#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SpatialModelVisitors

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/ip/diffim/BasisSets.h"
#include "lsst/ip/diffim/SpatialModelVisitors.h"

using lsst::pex::policy::Policy;
using lsst::pex::logging::Trace;

namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
using namespace lsst::ip::diffim;

typedef float PixelT;
typedef afwImage::Image<afwMath::Kernel::Pixel> ImageT;

BOOST_AUTO_TEST_CASE(KernelSumVisitor) {
    Policy::Ptr policy(new Policy);
    policy->set("kernelSumClipping", false);
    policy->set("maxKsumSigma", 3.0);

    /* Dummy images */
    afwImage::MaskedImage<PixelT>::Ptr mimg1(
        new afwImage::MaskedImage<PixelT>(5,5)
        );
    afwImage::MaskedImage<PixelT>::Ptr mimg2(
        new afwImage::MaskedImage<PixelT>(5,5)
        );

    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand1(new KernelCandidate<PixelT>(10., 10., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand2(new KernelCandidate<PixelT>(20., 20., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand3(new KernelCandidate<PixelT>(30., 30., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand4(new KernelCandidate<PixelT>(40., 40., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand5(new KernelCandidate<PixelT>(50., 50., mimg1, mimg2));
    
    ImageT img1(10,10);
    img1 = 0;
    *img1.at(4, 4) = 1;
    ImageT img2(10,10);
    img2 = 0;
    *img2.at(4, 4) = 1;
    ImageT img3(10,10);
    img3 = 0;
    *img3.at(4, 4) = 1;
    ImageT img4(10,10);
    img4 = 0;
    *img4.at(4, 4) = 1;
    /* outlier */
    ImageT img5(10,10);
    img5 = 0;
    *img5.at(4, 4) = 100;
    
    afwMath::Kernel::Ptr k1(new afwMath::FixedKernel(img1));
    afwMath::Kernel::Ptr k2(new afwMath::FixedKernel(img2));
    afwMath::Kernel::Ptr k3(new afwMath::FixedKernel(img3));
    afwMath::Kernel::Ptr k4(new afwMath::FixedKernel(img4));
    afwMath::Kernel::Ptr k5(new afwMath::FixedKernel(img5));

    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->setKernel(k1);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand2))->setKernel(k2);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand3))->setKernel(k3);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand4))->setKernel(k4);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand5))->setKernel(k5);

    /* Lets assume they got here being labeled as GOOD */
    cand1->setStatus(afwMath::SpatialCellCandidate::GOOD);
    cand2->setStatus(afwMath::SpatialCellCandidate::GOOD);
    cand3->setStatus(afwMath::SpatialCellCandidate::GOOD);
    cand4->setStatus(afwMath::SpatialCellCandidate::GOOD);
    cand5->setStatus(afwMath::SpatialCellCandidate::GOOD);
    
    /* Basic sum */
    {
        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand5));
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
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processCandidate(&(*cand3));
        kernelSumVisitor.processCandidate(&(*cand4));
        kernelSumVisitor.processCandidate(&(*cand5));
        kernelSumVisitor.processKsumDistribution();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processCandidate(&(*cand3));
        kernelSumVisitor.processCandidate(&(*cand4));
        kernelSumVisitor.processCandidate(&(*cand5));
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumMean(),  1.0);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumNpts(),  4);   /* cand2 sigma clipped out in the stats */
        BOOST_CHECK_EQUAL(kernelSumVisitor.getNRejected(), 0);   /* but not from the cells */
        BOOST_CHECK_EQUAL(cand1->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand2->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand3->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand4->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand5->getStatus(), afwMath::SpatialCellCandidate::GOOD);
    }

    /* There is sigma clipping in the stats, and now kernelSumClipping = true */
    {
        policy->set("kernelSumClipping", true);
        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processCandidate(&(*cand3));
        kernelSumVisitor.processCandidate(&(*cand4));
        kernelSumVisitor.processCandidate(&(*cand5));
        kernelSumVisitor.processKsumDistribution();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
        kernelSumVisitor.processCandidate(&(*cand1));
        kernelSumVisitor.processCandidate(&(*cand2));
        kernelSumVisitor.processCandidate(&(*cand3));
        kernelSumVisitor.processCandidate(&(*cand4));
        kernelSumVisitor.processCandidate(&(*cand5));
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumMean(),  1.0);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumNpts(),  4);   /* cand5 sigma clipped out in the stats */
        BOOST_CHECK_EQUAL(kernelSumVisitor.getNRejected(), 1);   /* and from the cells */
        BOOST_CHECK_EQUAL(cand1->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand2->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand3->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand4->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand5->getStatus(), afwMath::SpatialCellCandidate::BAD);
    }

    /* Finally test it in the visitCandidates() paradigm */
    {
        /* Reset its status */
        cand5->setStatus(afwMath::SpatialCellCandidate::UNKNOWN);

        afwMath::SpatialCellSet kernelCells = 
            afwMath::SpatialCellSet(afwImage::BBox(afwImage::PointI(0,0), 100, 100), 10, 10);
        kernelCells.insertCandidate(cand1);
        kernelCells.insertCandidate(cand2);
        kernelCells.insertCandidate(cand3);
        kernelCells.insertCandidate(cand4);
        kernelCells.insertCandidate(cand5);

        detail::KernelSumVisitor<PixelT> kernelSumVisitor(*policy);
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::AGGREGATE);
        kernelCells.visitCandidates(&kernelSumVisitor, 1);
        kernelSumVisitor.processKsumDistribution();
        kernelSumVisitor.setMode(detail::KernelSumVisitor<PixelT>::REJECT);
        kernelCells.visitCandidates(&kernelSumVisitor, 1);

        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumMean(),  1.0);   
        BOOST_CHECK_EQUAL(kernelSumVisitor.getkSumNpts(),  4);   /* cand5 sigma clipped out in the stats */
        BOOST_CHECK_EQUAL(kernelSumVisitor.getNRejected(), 1);   /* and from the cells */
        BOOST_CHECK_EQUAL(cand1->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand2->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand3->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand4->getStatus(), afwMath::SpatialCellCandidate::GOOD);
        BOOST_CHECK_EQUAL(cand5->getStatus(), afwMath::SpatialCellCandidate::BAD);
    }
}

BOOST_AUTO_TEST_CASE(SetPcaImageVisitor) {
    /* Dummy images */
    afwImage::MaskedImage<PixelT>::Ptr mimg1(
        new afwImage::MaskedImage<PixelT>(100,100)
        );
    afwImage::MaskedImage<PixelT>::Ptr mimg2(
        new afwImage::MaskedImage<PixelT>(100,100)
        );

    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand1(new KernelCandidate<PixelT>(10., 10., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand2(new KernelCandidate<PixelT>(20., 20., mimg1, mimg2));
    afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand3(new KernelCandidate<PixelT>(30., 30., mimg1, mimg2));
    
    ImageT img1(10,10);
    img1 = 0;
    *img1.at(4, 4) = 1.;
    ImageT img2(10,10);
    img2 = 0;
    *img2.at(4, 4) = 1.;
    /* Outlier */
    ImageT img3(10,10);
    img3 = 0;
    *img3.at(4, 4) = 100.;
    
    afwMath::Kernel::Ptr k1(new afwMath::FixedKernel(img1));
    afwMath::Kernel::Ptr k2(new afwMath::FixedKernel(img2));
    afwMath::Kernel::Ptr k3(new afwMath::FixedKernel(img3));

    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->setKernel(k1);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand2))->setKernel(k2);
    dynamic_cast<KernelCandidate<PixelT> *>(&(*cand3))->setKernel(k3);

    /* Test eigen decomposition; all images the same shape */
    {
        afwImage::ImagePca<ImageT> imagePca;
        detail::SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
        importStarVisitor.processCandidate(&(*cand1));
        importStarVisitor.processCandidate(&(*cand2));
        importStarVisitor.processCandidate(&(*cand3));

        imagePca.analyze();
        std::vector<ImageT::Ptr> eigenImages = imagePca.getEigenImages();
        std::vector<double> eigenValues = imagePca.getEigenValues();
        BOOST_CHECK_EQUAL(static_cast<int>(eigenImages.size()), 3);   
        BOOST_CHECK_CLOSE(eigenValues[0], 1., 1e-6);   
        BOOST_CHECK_SMALL(eigenValues[1], 1e-6);   
        BOOST_CHECK_SMALL(eigenValues[2], 1e-6);   
    }

    /* Test mean subtraction; recall all images are scaled to have the same kSum! */
    {
        afwImage::ImagePca<ImageT> imagePca;
        detail::SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
        importStarVisitor.processCandidate(&(*cand1));
        importStarVisitor.processCandidate(&(*cand2));
        importStarVisitor.processCandidate(&(*cand3));

        /* Have we divided by the kernel sum? */
        afwImage::ImagePca<ImageT>::ImageList imageList = imagePca.getImageList();
        BOOST_CHECK_CLOSE((*imageList[0])(4, 4), 1., 1e-6);
        BOOST_CHECK_CLOSE((*imageList[1])(4, 4), 1., 1e-6);
        BOOST_CHECK_CLOSE((*imageList[2])(4, 4), 1., 1e-6);

        importStarVisitor.subtractMean();
        ImageT kMean = *(importStarVisitor.returnMean());

        /* Mean has the right central pixel value */
        BOOST_CHECK_CLOSE(kMean(4, 4), 1., 1e-6);

        /* Candidates in imagePca have their values modified before doing Pca */
        imageList = imagePca.getImageList();
        BOOST_CHECK_SMALL((*imageList[0])(4, 4), 1e-6);
        BOOST_CHECK_SMALL((*imageList[1])(4, 4), 1e-6);
        BOOST_CHECK_SMALL((*imageList[2])(4, 4), 1e-6);

    }
}

BOOST_AUTO_TEST_CASE(BuildSingleKernelVisitor) {
    Policy::Ptr policy(new Policy);
    policy->set("constantVarianceWeighting", false);
    policy->set("iterateSingleKernel", false);
    policy->set("singleKernelClipping", true);
    policy->set("candidateResidualMeanMax", 0.25);
    policy->set("candidateResidualStdMax", 1.25);

    int imageSize = 35;
    int kernelSize = 11;
    lsst::afw::math::KernelList basisList = generateDeltaFunctionBasisSet(kernelSize, kernelSize);
    PsfMatchingFunctor<PixelT> kFunctor = PsfMatchingFunctor<PixelT>(basisList);

    afwMath::GaussianFunction2<double> gaussFunc2(2, 2);
    afwMath::GaussianFunction2<double> gaussFunc4(4, 4);
    afwMath::AnalyticKernel k1(imageSize, imageSize, gaussFunc2);
    afwMath::AnalyticKernel k2(imageSize, imageSize, gaussFunc4);
    afwImage::Image<double> kImage1D(k1.getDimensions());
    afwImage::Image<double> kImage2D(k2.getDimensions());
    (void)k1.computeImage(kImage1D, true);
    (void)k2.computeImage(kImage2D, true);

    {
        /* Make them the same so a delta function result */
        afwImage::Image<PixelT> kImage1F(kImage1D, true);
        afwImage::Image<PixelT> kImage2F(kImage1D, true);
        afwImage::MaskedImage<PixelT>::Ptr mimg1(
            new afwImage::MaskedImage<PixelT>(k1.getDimensions())
            );
        afwImage::MaskedImage<PixelT>::Ptr mimg2(
            new afwImage::MaskedImage<PixelT>(k2.getDimensions())
            );
        *mimg1->getImage() = kImage1F;
        *mimg2->getImage() = kImage2F;
        *mimg1->getVariance() = 1;
        *mimg2->getVariance() = 1;
        *mimg1->getMask() = 0x0;
        *mimg2->getMask() = 0x0;

        afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand1(new KernelCandidate<PixelT>(10., 10., mimg1, mimg2));
        BOOST_CHECK(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->hasKernel() == false);

        detail::BuildSingleKernelVisitor<PixelT> singleKernelFitter(kFunctor, *policy);
        singleKernelFitter.processCandidate(&(*cand1));

        BOOST_CHECK(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->hasKernel() == true);
        /* We're not checking kFunctor here, just that the values are set, so just check to 1% */
        BOOST_CHECK_CLOSE(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getKsum(), 1.0, 1.);
        BOOST_CHECK_SMALL(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getBackground(), 1.);
        BOOST_CHECK_EQUAL(cand1->getStatus(), afwMath::SpatialCellCandidate::GOOD);
    }

    /* Check to make sure we dont overwrite kernel and background,
     * but we do overwrite M and B */
    {
        /* Make them the same so a delta function result */
        afwImage::Image<PixelT> kImage1F(kImage1D, true);
        afwImage::Image<PixelT> kImage2F(kImage1D, true);
        afwImage::MaskedImage<PixelT>::Ptr mimg1(
            new afwImage::MaskedImage<PixelT>(k1.getDimensions())
            );
        afwImage::MaskedImage<PixelT>::Ptr mimg2(
            new afwImage::MaskedImage<PixelT>(k2.getDimensions())
            );
        *mimg1->getImage() = kImage1F;
        *mimg2->getImage() = kImage2F;
        *mimg1->getVariance() = 1;
        *mimg2->getVariance() = 1;
        *mimg1->getMask() = 0x0;
        *mimg2->getMask() = 0x0;
        afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand1(new KernelCandidate<PixelT>(10., 10., mimg1, mimg2));
        BOOST_CHECK(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->hasKernel() == false);
        detail::BuildSingleKernelVisitor<PixelT> singleKernelFitter(kFunctor, *policy);
        singleKernelFitter.processCandidate(&(*cand1));

        /* Reprocess with a new basis */
        Eigen::MatrixXd Morig = *(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getM());
        Eigen::VectorXd Borig = *(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getB());
        lsst::afw::math::Kernel::Ptr Korig = dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getKernel();
        double bgorig                      = dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getBackground();

        lsst::afw::math::KernelList basisList2 = generateDeltaFunctionBasisSet(kernelSize-1, kernelSize-1);
        PsfMatchingFunctor<PixelT> kFunctor2 = PsfMatchingFunctor<PixelT>(basisList2);
        detail::BuildSingleKernelVisitor<PixelT> singleKernelFitter2(kFunctor2, *policy);
        singleKernelFitter2.setCandidateKernel(false);
        singleKernelFitter2.setSkipBuilt(false);
        singleKernelFitter2.processCandidate(&(*cand1));
        
        Eigen::MatrixXd Mnew = *(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getM());
        Eigen::VectorXd Bnew = *(dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getB());
        lsst::afw::math::Kernel::Ptr Knew  = dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getKernel();
        double bgnew                       = dynamic_cast<KernelCandidate<PixelT> *>(&(*cand1))->getBackground();

        /* Matrices are updated for spatial fitting */
        BOOST_CHECK(Borig.size() == (kernelSize*kernelSize+1));
        BOOST_CHECK(Bnew.size()  == ((kernelSize-1)*(kernelSize-1)+1));
        BOOST_CHECK(Morig.rows() > Mnew.rows());
        BOOST_CHECK(Morig.cols() > Mnew.cols());
        /* But the raw kernel is not */
        BOOST_CHECK(Korig->getHeight() == Knew->getHeight());
        BOOST_CHECK(Korig->getWidth() == Knew->getWidth());
        BOOST_CHECK(bgorig == bgnew);
    }
    
}

BOOST_AUTO_TEST_CASE(BuildSpatialKernelVisitor) {
    int spatialKernelOrder = 2;
    int spatialBgOrder = 1;

    Policy::Ptr policy(new Policy);
    policy->set("spatialKernelOrder", spatialKernelOrder);
    policy->set("spatialBgOrder", spatialBgOrder);
    policy->set("kernelBasisSet", "delta-function");
    policy->set("usePcaForSpatialKernel", true);


    int kernelSize = 11;
    lsst::afw::math::KernelList basisList = generateDeltaFunctionBasisSet(kernelSize, kernelSize);
    detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(basisList, spatialKernelOrder, spatialBgOrder, *policy);

    {
        /* This is going to be tough to unit test thoroughly */
    }
}

BOOST_AUTO_TEST_CASE(AssessSpatialKernelVisitor) {
    int spatialKernelOrder = 2;
    int spatialBgOrder = 1;
    unsigned int kSize = 5;
    unsigned int nBases = kSize * kSize;

    Policy::Ptr policy(new Policy);
    policy->set("spatialKernelOrder", spatialKernelOrder);
    policy->set("spatialBgOrder", spatialBgOrder);
    policy->set("spatialKernelClipping", true);
    policy->set("candidateResidualMeanMax", 0.25);
    policy->set("candidateResidualStdMax", 1.25);

    lsst::afw::math::KernelList basisList = generateDeltaFunctionBasisSet(kSize, kSize);
    /* Spatial Kernel */
    afwMath::Kernel::SpatialFunctionPtr spatialKernelFunction(
        new afwMath::PolynomialFunction2<double>(spatialKernelOrder)
        );
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    for (unsigned int i = 0; i < nBases; i++) {
        afwMath::Kernel::SpatialFunctionPtr spatialFunction(spatialKernelFunction->copy());
        spatialFunctionList.push_back(spatialFunction);
    }
    afwMath::LinearCombinationKernel::Ptr spatialKernel(new afwMath::LinearCombinationKernel(basisList, spatialFunctionList));
    
    /* Spatial Background */
    afwMath::Kernel::SpatialFunctionPtr spatialBg(
        new afwMath::PolynomialFunction2<double>(spatialBgOrder)
        );

    /* Set up some fake terms */
    std::vector<std::vector<double> > kCoeffs;
    kCoeffs.reserve(nBases);
    for (unsigned int i = 0, idx = 0; i < nBases; i++) {
        kCoeffs.push_back(std::vector<double>(spatialKernelFunction->getNParameters()));
        for (unsigned int j = 0; j < spatialKernelFunction->getNParameters(); j++) {
            kCoeffs[i][j] = sqrt(idx++);
        }
    }
    spatialKernel->setSpatialParameters(kCoeffs);

    /* Set up fake background coefficients; all zero */
    std::vector<double> bgCoeffs(spatialBg->getNParameters());
    spatialBg->setParameters(bgCoeffs);

    {
        Trace::setVerbosity("lsst.ip.diffim", 6);
        /* Manually create a canidate with the spatial kernel, and
         * then assess its quality with the same spatial kernel */

        unsigned int loc = 50;
        afwImage::MaskedImage<PixelT>::Ptr mimg1(
            new afwImage::MaskedImage<PixelT>(100,100)
            );
        *mimg1->at(loc, loc) = afwImage::MaskedImage<PixelT>::Pixel(1, 0x0, 1);
        afwImage::MaskedImage<PixelT>::Ptr mimg2(
            new afwImage::MaskedImage<PixelT>(mimg1->getDimensions())
            );
        afwMath::convolve(*mimg2, *mimg1, *spatialKernel, false);

        /* Grab a subimage to create the candidate */
        afwImage::BBox bbox = afwImage::BBox(afwImage::PointI(loc-10, loc-10),
                                             afwImage::PointI(loc+10, loc+10));
        afwImage::MaskedImage<PixelT>::Ptr tmi(
            new afwImage::MaskedImage<PixelT>(*mimg1, bbox)
            );
        afwImage::MaskedImage<PixelT>::Ptr smi(
            new afwImage::MaskedImage<PixelT>(*mimg2, bbox)
            );

        tmi->writeFits("/tmp/tmi");
        smi->writeFits("/tmp/smi");

        /* Create the candidate and give it the spatial kernel
         * evaluated at its position */
        afwMath::SpatialCellImageCandidate<ImageT>::Ptr cand(new KernelCandidate<PixelT>(loc, loc, tmi, smi));
        afwImage::Image<double> kImage(spatialKernel->getDimensions());
        (void)spatialKernel->computeImage(kImage, false, loc, loc);
        boost::shared_ptr<afwMath::Kernel>
            kernelPtr(new afwMath::FixedKernel(kImage));
        kImage.writeFits("/tmp/kernel1.fits");
        dynamic_cast<KernelCandidate<PixelT> *>(&(*cand))->setKernel(kernelPtr);

        /* Process it! */
        detail::AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, spatialBg, *policy);
        spatialKernelAssessor.processCandidate(&(*cand));
    }
}

