#include <lsst/pex/logging/Trace.h>
#include <lsst/daf/base/Citizen.h>
#include <lsst/afw/MaskedImage.h>
#include <lsst/detection/Footprint.h>

using namespace std;
using namespace lsst::afw;

typedef lsst::afw::maskPixelType MaskT;
typedef float ImageT;
typedef double KernelT;
typedef double FuncT;

int main( int argc, char** argv )
{
    {
        lsst::pex::logging::Trace::setDestination(cout);
        lsst::pex::logging::Trace::setVerbosity(".", 4);

        string templateImage = argv[1];
        MaskedImage<ImageT,MaskT> templateMaskedImage;
        try {
            templateMaskedImage.readFits(templateImage);
        } catch (lsst::pex::exceptions::ExceptionStack &e) {
            cerr << "Failed to open template image " << templateImage << ": " << e.what() << endl;
            return 1;
        }

        float threshold = atof(argv[2]);

        // Find detections
        lsst::detection::DetectionSet<ImageT,MaskT> 
            detectionSet(templateMaskedImage, lsst::detection::Threshold(threshold, lsst::detection::Threshold::VALUE));
        vector<lsst::detection::Footprint::PtrType> footprintVector = detectionSet.getFootprints();
        cout << " Detected " << footprintVector.size() << " footprints at value threshold " << threshold << endl;

    }
    
    if (Citizen::census(0) == 0) {
        cerr << "No leaks detected" << endl;
    } else {
        cerr << "Leaked memory blocks:" << endl;
        Citizen::census(cerr);
    } 
}
