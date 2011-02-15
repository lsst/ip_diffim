#include "Eigen/Core"
#include "lsst/afw/math.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/ip/diffim.h"

namespace ipDiffim       = lsst::ip::diffim;
namespace afwMath        = lsst::afw::math;
namespace afwGeom        = lsst::afw::geom;
namespace afwImage       = lsst::afw::image;
int main() {
    int dimen = 11;
    afwImage::Image<float> foo(dimen,dimen);
    for (int y = 0, n = 0; y < dimen; y++) {
        for (int x = 0; x < dimen; x++, n++) {
            foo(x,y) = n;
        }
    }
    foo.writeFits("maskedKernel.fits");

    Eigen::MatrixXd test = ipDiffim::imageToEigenMatrix(foo);
    /*
    test << 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
        80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 
        70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
        60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
        50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 
        40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
        30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
    */

    afwImage::Image<double> kImage(5, 5);
    afwMath::FixedKernel kernel(kImage);

    std::cout << "All data:" << std::endl << test << std::endl << std::endl;
    int startCol = kernel.getCtrX();
    int startRow = kernel.getCtrY();
    int endCol = dimen - (kernel.getWidth() - kernel.getCtrX());
    int endRow = dimen - (kernel.getHeight() - kernel.getCtrY());
    endCol += 1;
    endRow += 1;
    
    std::cout << "A " << startCol << " " << startRow << " " << endCol << " " << endRow << std::endl;
        
    Eigen::MatrixXd subimage = test.block(startRow, 
                                          startCol, 
                                          endRow-startRow, 
                                          endCol-startCol);
    std::cout << "Good pixels after convolution:" << std::endl << subimage << std::endl << std::endl;

    /* We want to ignore all pixels from 4,4 to 6,6 */
    afwGeom::BoxI maskBox(afwGeom::makePointI(4, 4),
                          afwGeom::makePointI(6, 6));
    
    int maskStartCol = maskBox.getMinX();
    int maskStartRow = maskBox.getMinY();
    int maskEndCol   = maskBox.getMaxX();
    int maskEndRow   = maskBox.getMaxY();
    std::cout << "B " << maskBox.getMinX() << " " << maskBox.getMinY() << " " << maskBox.getHeight() << " " << maskBox.getWidth() << std::endl;
    std::cout << "B " << maskStartCol << " " << maskStartRow << " " << maskEndCol << " " << maskEndRow << std::endl;

    /*
    afwGeom::BoxI tBox = afwGeom::BoxI(afwGeom::makePointI(startCol, maskEndRow),
                                       afwGeom::makePointI(endCol, endRow));
    
    afwGeom::BoxI bBox = afwGeom::BoxI(afwGeom::makePointI(startCol, startRow),
                                       afwGeom::makePointI(endCol, maskStartRow));
    
    afwGeom::BoxI lBox = afwGeom::BoxI(afwGeom::makePointI(startCol, maskStartRow),
                                       afwGeom::makePointI(maskStartCol, maskEndRow));
    
    afwGeom::BoxI rBox = afwGeom::BoxI(afwGeom::makePointI(maskEndCol, maskStartRow),
                                       afwGeom::makePointI(endCol, maskEndRow));
    */

    endCol -= 1;
    endRow -= 1;
    afwImage::BBox tBox = afwImage::BBox(afwImage::PointI(startCol, maskEndRow + 1),
                                         afwImage::PointI(endCol, endRow));
    
    afwImage::BBox bBox = afwImage::BBox(afwImage::PointI(startCol, startRow),
                                         afwImage::PointI(endCol, maskStartRow - 1));
    
    afwImage::BBox lBox = afwImage::BBox(afwImage::PointI(startCol, maskStartRow),
                                         afwImage::PointI(maskStartCol - 1, maskEndRow));
    
    afwImage::BBox rBox = afwImage::BBox(afwImage::PointI(maskEndCol + 1, maskStartRow),
                                         afwImage::PointI(endCol, maskEndRow));


    std::vector<afwImage::BBox> boxArray;
    boxArray.push_back(tBox);
    boxArray.push_back(bBox);
    boxArray.push_back(lBox);
    boxArray.push_back(rBox);
    std::vector<afwImage::BBox>::iterator biter = boxArray.begin();
    for (; biter != boxArray.end(); ++biter) {
        int area = (*biter).getWidth() * (*biter).getHeight();
        afwImage::Image<float> subimage = afwImage::Image<float>(foo, (*biter));
        Eigen::MatrixXd subeigen = ipDiffim::imageToEigenMatrix(subimage);
        std::cout << area << std::endl;
        std::cout << subeigen << std::endl << std::endl;
        subeigen.resize(area, 1);
        std::cout << subeigen << std::endl << std::endl;
    }
}

        
    
