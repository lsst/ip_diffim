/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
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
    int vals111r0[] = {-4, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; vals111[0] = vals111r0;
    int vals111r1[] = {1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}; vals111[1] = vals111r1;
    int vals111r2[] = {0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0}; vals111[2] = vals111r2;
    int vals111r3[] = {1, 0, 1, -4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}; vals111[3] = vals111r3;
    int vals111r4[] = {1, 0, 0, 0, -4, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0}; vals111[4] = vals111r4;
    int vals111r5[] = {0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0}; vals111[5] = vals111r5;
    int vals111r6[] = {0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0}; vals111[6] = vals111r6;
    int vals111r7[] = {0, 0, 0, 1, 1, 0, 1, -4, 0, 0, 0, 1, 0, 0, 0, 0}; vals111[7] = vals111r7;
    int vals111r8[] = {0, 0, 0, 0, 1, 0, 0, 0, -4, 1, 0, 1, 1, 0, 0, 0}; vals111[8] = vals111r8;
    int vals111r9[] = {0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0}; vals111[9] = vals111r9;
    int vals111r10[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0}; vals111[10] = vals111r10;
    int vals111r11[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, -4, 0, 0, 0, 1}; vals111[11] = vals111r11;
    int vals111r12[] = {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -4, 1, 0, 1}; vals111[12] = vals111r12;
    int vals111r13[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0}; vals111[13] = vals111r13;
    int vals111r14[] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1}; vals111[14] = vals111r14;
    int vals111r15[] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, -4}; vals111[15] = vals111r15;
    knownBList.push_back(vals111);

    // - order 1, tapered, central
    Param pars121 = {1, ipDiffim::TAPERED, ipDiffim::CENTRAL_DIFFERENCE};
    paramSetList.push_back(pars121);
    int *vals121[16];
    int vals121r0[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[0] = vals121r0;
    int vals121r1[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[1] = vals121r1;
    int vals121r2[] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[2] = vals121r2;
    int vals121r3[] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[3] = vals121r3;
    int vals121r4[] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[4] = vals121r4;
    int vals121r5[] = {0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0}; vals121[5] = vals121r5;
    int vals121r6[] = {0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0, 0, 0, 0}; vals121[6] = vals121r6;
    int vals121r7[] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}; vals121[7] = vals121r7;
    int vals121r8[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}; vals121[8] = vals121r8;
    int vals121r9[] = {0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0, 0}; vals121[9] = vals121r9;
    int vals121r10[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -4, 1, 0, 0, 1, 0}; vals121[10] = vals121r10;
    int vals121r11[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}; vals121[11] = vals121r11;
    int vals121r12[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}; vals121[12] = vals121r12;
    int vals121r13[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}; vals121[13] = vals121r13;
    int vals121r14[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}; vals121[14] = vals121r14;
    int vals121r15[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; vals121[15] = vals121r15;
    knownBList.push_back(vals121);
    

    // Loop over each possible set of parameters and check the outputs
    for (unsigned int k = 0; k < paramSetList.size(); k++) {
        
        int const order                 = paramSetList[k].order;
        ipDiffim::BoundStyle boundStyle = paramSetList[k].bound;
        ipDiffim::DiffStyle  diffStyle  = paramSetList[k].diff;
        int **bKnown = knownBList[k];

        // note: calling FdrBMatrix (which returns B),
        //       not FiniteDifferenceRegularization (which returns h = BTranspose*B)
        boost::shared_ptr<Eigen::MatrixXd> b =
            ipDiffim::details::generateFdrBMatrix(nCols, nRows, order, boundStyle, diffStyle);

        // Do the comparison
        for (int i = 0; i < nRows*nCols; i++) {
            for (int j = 0; j < nRows*nCols; j++) {
                int val = (*b)(i, j);
                int valKnown = bKnown[i][j];
                BOOST_CHECK_EQUAL(val, valKnown);
            }
        }

        
#if defined(PRINT2D)
        // print the matrix elements in 2D, so you can see how they're positioned
        // -- use this output to verify by eye
        for (int i = 0; i < nRows*nCols; i++) {
            for (int j = 0; j < nRows*nCols; j++) {
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
        for (int i = 0; i < nRows*nCols; i++) {
            std::cout << "    int vals" << order << boundStyle << diffStyle << "_" << i << "[] = {";
            for (int j = 0; j < nRows*nCols; j++) {
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
