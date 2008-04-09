// -*- lsst-c++ -*-
#include "lsst/afw/Trace.h"
#include "lsst/afw/math/Kernel.h"

using namespace std;
using namespace lsst::afw;

int main( int argc, char** argv )
{
    // ENTIRELY PSEUDO-CODE

    Trace::setDestination(cout);
    Trace::setVerbosity(".", 0);
    
    typedef uint8 MaskPixelType;
    typedef float32 ImagePixelType;

    // This will create the delta function kernel collection; do once
    vector<Kernel> deltaFunctionSet;
    int sizeKernel = 8;
    for (unsigned row = 0; row < sizeKernel; row++) {
        for (unsigned col = 0; col < sizeKernel; col++) {
            // This does not exist yet
            DeltaFunctionKernel<ImagePixelType> kernel(sizeKernel, sizeKernel, row, col);
            Image<ImagePixelType> kImage(kernel.getImage());
            kImage.writeFits();
            deltaFunctionSet.push_back(kernel);
        }
    }

    // This will get the coefficient for each basis kernel; do for each object
    // Need object collection here!
    // Collection<Object> objectCollection();
    vector<Object> objectCollection;

    // Reusable view around each object
    MaskedImage<ImagePixelType,MaskPixelType>::MaskedImagePtrT convolvePtr;
    MaskedImage<ImagePixelType,MaskPixelType>::MaskedImagePtrT nonconvolvePtr;

    // Should be a Collection; use STL for now
    vector<Image> basisImageCollection;

    // Iterate over object; use iterator instead?
    for (unsigned nobj = 0; nobj < objectCollection.size(); nobj++) {
        Object<> diffImObject = objectCollection[nobj];

        // grab view around each object
        // do i really want a new stamp or just a view?
        BBox2i stamp(diffImObject.rowc - diffImObject.drow, diffImObject.rowc + diffImObject.drow, 
                     diffImObject.colc - diffImObject.dcol, diffImObject.colc + diffImObject.dcol);

        convolvePtr    = convolveMaskedImage.getSubImage(stamp);
        nonconvolvePtr = nonconvolveMaskedImage.getSubImage(stamp);
        
        vector B;
        matrix M;

        // integral over image's dx and dy
        for (unsigned row = 0; row < convolvePtr->getImage().nrow; row++) {
            for (unsigned col = 0; col < convolvePtr->getImage().ncol; col++) {
                
                ncData = nonconvolvePtr[row][col];
                
                // kernel index i
                for (unsigned kerni = 0; kerni < deltaFunctionSet.size(); kerni++) {
                    kerneli = deltaFunctioSet[kerni];
                    cData   = convolvePtr[row][col];
                    // This is effectively what we want to happen with the delta function kernel.
                    // cData = convolvePtr[row + kerneli->drow][col + kerneli->dcol];
                    // The more general case is

                    // convolve
                    boost::shared_ptr<MaskedImage<ImagePixelType, MaskPixelType> >
                        cImagePtri = Kernel::convolve(convolvePtr[row][col], kerneli, threshold, vw::NoEdgeExtension(), -1);
            
                    B[kerneli] += ncData.camera() * cImagePtri.camera() / (ncData.variance() + cImagePtri.variance());

                    // background term (is this cData right?)
                    M[kerneli][nkernel] += cImagePtri.camera() / (cData.variance() + cImagePtri.variance());
                    
                    // kernel index j 
                    for (unsigned kernj = 0; kernj < deltaFunctionSet.size(); kernj++) {
                        kernelj = deltaFunctionSet[kernj];
                        
                        // convolve
                        boost::shared_ptr<MaskedImage<ImagePixelType, MaskPixelType> >
                            cImagePtrj = Kernel::convolve(convolvePtr[row][col], kernelj, threshold, vw::NoEdgeExtension(), -1);
                    
                        M[kerni][kernj] += cImagePtri.camera() * cImagePtrj.camera() / (cImagePtri.variance() * cImagePtrj.variance());
                    }
                }
                
                // background
                M[nkernel][nkernel] += 1. / (ncData.variance() + cData.variance());
                B[knernel]          += ncData.camera() / (ncData.variance() + cData.variance());
            } // col
        } // row
        
        vector LUv;
        matrix LUm = LUD(LUv, M);

        vector solution;
        LUSolve(LUm, B, LUv);

        // rearrange vector into kernel Image
        image KernelImage;

        basisImageCollection.push_back(KernelImage);
    }
    
    nreplace = 1;
    niter = 0;
    while ( (niter < policy.maxIterPCA) && (nreplace > 0) ) {
        nreplace = 0;

        // Run SVD on basisImageCollection to get principal components
        vector<Image> principalComponents;
        vector<ImagePixelType> principalCoefficients (basisImageCollection.size(), policy.nPrincipalComponents);

        int nComponents = 3;
        for (unsigned nbasis = 0; nbasis < basisImageCollection.size(); nbasis++) {
            MaskedImage Reconstruction;
            for (unsigned component = 0; component < nComponents; component++) {
                principalCoefficients[nbasis, component] = basisImageCollection[nbasis] * principalComponents[component];
                Reconstruction += principalCoefficients[nbasis, component] * basisImageCollection[nbasis];
            }
            MaskedImage difference = basisImageCollection[nbasis] - Reconstruction;
            
            float chi2 = difference.chi2();
            
            if (chi2 > policy.chi2PCA) {
                // remove from basisImageCollection;
                nreplace += 1;
            }
        }
    }
    // now we have the good set of principalCoefficients + principalComponents
    // fit to spatially varying function

    // and you are done.
}
