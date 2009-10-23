// -*- lsst-c++ -*-

#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SpatialModelQuality

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

typedef float Pixel;
typedef afwImage::Image<afwMath::Kernel::Pixel> Image;

//Trace::setVerbosity("lsst.ip.diffim", 6);

void JackknifeResample(afwMath::SpatialCellSet kernelCells) {
    std::vector<int> goodIds;
    afwMath::CellList cellList = kernelCells.getCellList();
    for (afwMath::CellList::iterator cell = cellList.begin(), end = cellList.end(); cell != end; ++cell) {
        for (SpatialCell::iterator candidate = (*cell)->begin(), candidateEnd = (*cell)->end();
             candidate != candidateEnd; ++candidate) {
            if ((*candidate)->getStatus() == afwMath::SpatialCell::GOOD) {
                goodIds.push_back((*candidate)->getId());
            }
        }
    }
    int nCandidates = goodIds.size();

    policy.set("spatialKernelClipping", false);
    for (int i = 0; i < nCandidates; i++) {
        /* Remove this object from the fitting */
        kernelCells.getCandidateById(i).setStatus(afwMath::SpatialCell::BAD);
        double chi2a = kernelCells.getCandidateById(i).getChi2();

        /* Revisit the other stamps to rebuild the kernel */
        detail::BuildSpatialKernelVisitor<PixelT> spatialKernelFitter(*basisListToUse, 
                                                                      spatialKernelOrder, 
                                                                      spatialBgOrder, 
                                                                      policy);
        kernelCells.visitCandidates(&spatialKernelFitter, nStarPerCell);
        spatialKernelFitter.solveLinearEquation();
        std::pair<afwMath::LinearCombinationKernel::Ptr, 
            afwMath::Kernel::SpatialFunctionPtr> KB = spatialKernelFitter.getSpatialModel();
        spatialKernel     = KB.first;
        spatialBackground = KB.second;
        nFit              = spatialKernelFitter.getNCandidates();
        BOOST_CHECK_EQUAL(nFit, nCandidates-1);

        /* Assess the quality of the one stamp we removed */
        detail::AssessSpatialKernelVisitor<PixelT> spatialKernelAssessor(spatialKernel, 
                                                                         spatialBackground, 
                                                                         policy);
        spatialKernelAssessor.visitCandidate(kernelCells.getCandidateById(i));
        double chi2b = kernelCells.getCandidateById(i).getChi2();

        BOOST_CHECK((chi2b - chi2a) < 0.25);
        BOOST_CHECK(chi2b < 1.5);
        
        /* And, put it back for the next fit */
        kernelCells.getCandidateById(i).setStatus(afwMath::SpatialCell::GOOD);
    }
}

BOOST_AUTO_TEST_CASE(testDeltaFunction) {
}

BOOST_AUTO_TEST_CASE(testRegularizedDeltaFunction) {
}

BOOST_AUTO_TEST_CASE(testAlardLupton) {
}
