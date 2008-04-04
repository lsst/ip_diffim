#include <lsst/detection/Footprint.h>

using namespace std;
using namespace lsst::afw;

int main( int argc, char** argv )
{
    {
        double radius = 20;
        
        const vw::BBox2i region1(static_cast<int>(78.654-radius/2),
                                 static_cast<int>(3573.945-radius/2),
                                 static_cast<int>(radius), 
                                 static_cast<int>(radius));
        lsst::detection::Footprint::PtrType fp1(
            new lsst::detection::Footprint(1, region1)
            );
        
        BBox2i footprintBBox = fp1->getBBox();
        cout << "Reg1a " << region1.min() << " " << region1.max() << endl;
        cout << "Reg1b " << footprintBBox.min() << " " << footprintBBox.max() << endl;
    }
}
        
