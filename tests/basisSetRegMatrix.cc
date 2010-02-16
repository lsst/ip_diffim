#include <iostream>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasisSetRegMatrix

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "Eigen/Core"
#include "lsst/ip/diffim/BasisSets.h"

namespace ipDiffim = lsst::ip::diffim;

// #define PRINT2D 1
// #define PRINTCODE 1

typedef struct {
    int order;
    ipDiffim::BoundStyle bound;
    ipDiffim::DiffStyle diff;
} Param;

BOOST_AUTO_TEST_CASE(BasisSetRegMatrix) { /* parasoft-suppress  LsstDm-3-2a LsstDm-3-4a LsstDm-4-6 LsstDm-5-25 "Boost non-Std" */

    int nRows = 4, nCols = 4;

    std::vector<Param> paramSetList;     // stored as order, boundStyle, diffStyle
    std::vector<int**> knownBList;

    // check a few of the different regularization matrices against known values
    // These two have been checked by eye and are known to produce the correct regularization matrix
    // The issue is to verify that the tapering (lower order near kernel edge), and wrapping
    //    are correct.  The only way to do this is to print them as 2d matrices and look.

    // There's code below to print the matrices, and also to generate the code lines defining
    //    the test arrays, which have indeed been verified.
    // This test actually does relatively little ... if someone changes the regularization matrices
    //   (which is unlikely), it will raise a red flag.
    
    // boundStyles: 0 = unwrapped, 1 = wrapped, 2 = order-tapered ('order' is highest used)
    // diffStyles: 0 = forward, 1 = central

    
    // - order 1, wrapped, central
    Param pars111 = {1, ipDiffim::WRAPPED, ipDiffim::CENTRAL_DIFFERENCE};
    paramSetList.push_back(pars111);
    int *vals111[16];
    int vals111_0[] = {-4, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; vals111[0] = vals111_0;
    int vals111_1[] = {1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}; vals111[1] = vals111_1;
    int vals111_2[] = {0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0}; vals111[2] = vals111_2;
    int vals111_3[] = {1, 0, 1, -4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}; vals111[3] = vals111_3;
    int vals111_4[] = {1, 0, 0, 0, -4, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}; vals111[4] = vals111_4;
    int vals111_5[] = {0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0}; vals111[5] = vals111_5;
    int vals111_6[] = {0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0}; vals111[6] = vals111_6;
    int vals111_7[] = {0, 0, 0, 1, 1, 0, 1, -4, 0, 0, 0, 1, 0, 0, 0, 0}; vals111[7] = vals111_7;
    int vals111_8[] = {0, 0, 0, 0, 1, 0, 0, 0, -4, 1, 0, 1, 1, 0, 0, 0}; vals111[8] = vals111_8;
    int vals111_9[] = {0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0}; vals111[9] = vals111_9;
    int vals111_10[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0}; vals111[10] = vals111_10;
    int vals111_11[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, -4, 0, 0, 0, 1}; vals111[11] = vals111_11;
    int vals111_12[] = {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -4, 1, 0, 1}; vals111[12] = vals111_12;
    int vals111_13[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0}; vals111[13] = vals111_13;
    int vals111_14[] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1}; vals111[14] = vals111_14;
    int vals111_15[] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, -4}; vals111[15] = vals111_15;
    knownBList.push_back(vals111);

    // - order 1, tapered, central
    Param pars121 = {1, ipDiffim::TAPERED, ipDiffim::CENTRAL_DIFFERENCE};
    paramSetList.push_back(pars121);
    int *vals121[16];
    int vals121_0[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[0] = vals121_0;
    int vals121_1[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[1] = vals121_1;
    int vals121_2[] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[2] = vals121_2;
    int vals121_3[] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[3] = vals121_3;
    int vals121_4[] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[4] = vals121_4;
    int vals121_5[] = {0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0}; vals121[5] = vals121_5;
    int vals121_6[] = {0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0}; vals121[6] = vals121_6;
    int vals121_7[] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[7] = vals121_7;
    int vals121_8[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}; vals121[8] = vals121_8;
    int vals121_9[] = {0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0}; vals121[9] = vals121_9;
    int vals121_10[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0}; vals121[10] = vals121_10;
    int vals121_11[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}; vals121[11] = vals121_11;
    int vals121_12[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; vals121[12] = vals121_12;
    int vals121_13[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}; vals121[13] = vals121_13;
    int vals121_14[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}; vals121[14] = vals121_14;
    int vals121_15[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; vals121[15] = vals121_15;
    knownBList.push_back(vals121);
    

    // Loop over each possible set of parameters and check the outputs
    for(unsigned int k = 0; k < paramSetList.size(); k++) {
        
        int const order                 = paramSetList[k].order;
        ipDiffim::BoundStyle boundStyle = paramSetList[k].bound;
        ipDiffim::DiffStyle  diffStyle  = paramSetList[k].diff;
        int **bKnown = knownBList[k];
        boost::shared_ptr<Eigen::MatrixXd> b =
            ipDiffim::details::generateFdrBMatrix(nCols, nRows, order, boundStyle, diffStyle);

        // Do the comparison
        for(int i = 0; i < nRows*nCols; i++) {
            for(int j = 0; j < nRows*nCols; j++) {
                int val = (*b)(i, j);
                int valKnown = bKnown[i][j];
                BOOST_CHECK_EQUAL(val, valKnown);
            }
        }

        
#if defined(PRINT2D)
        // print the matrix elements in 2D, so you can see how they're positioned
        // -- use this output to verify by eye
        for(int i = 0; i < nRows*nCols; i++) {
            for(int j = 0; j < nRows*nCols; j++) {
                int val = (*b)(i, j);
                std::cout << val << " ";
                if (j % nCols == (nCols - 1)) std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
#endif

        
#if defined(PRINTCODE)
        // print out the function to return a hard-coded matrix
        // -- assuming you've verified the values by eye!!!
        std::cout << "    int *vals" << order << boundStyle << diffStyle <<
            "[" << nRows*nCols << "];" << std::endl;
        for(int i = 0; i < nRows*nCols; i++) {
            std::cout << "    int vals" << order << boundStyle << diffStyle << "_" << i << "[] = {";
            for(int j = 0; j < nRows*nCols; j++) {
                int val = (*b)(i, j);
                std::cout << val;
                if (j < nRows*nCols-1) std::cout << ", ";
            }
            std::cout << "}; vals" << order << boundStyle << diffStyle << "[" << i <<
                "] = vals" << order << boundStyle << diffStyle << "_" << i << ";"  << std::endl;
        }
#endif        

    }

}
