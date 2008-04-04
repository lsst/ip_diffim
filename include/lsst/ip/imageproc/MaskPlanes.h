// -*- lsst-c++ -*-
#ifndef LSST_Ip_Diffim_MaskPlanes_H
#define LSST_Ip_Diffim_MaskPlanes_H

#include <map>
#include <string>

typedef std::map<std::string, int> strIntMap;
typedef std::pair<std::string, int> strIntPair;

namespace lsst {
namespace ip {
namespace diffim {
        
    strIntMap MaskPlanes;
    MaskPlanes.insert(strIntPair{"BAD",   0}); /* For non-LSST data, not entirely sure why data are bad */
    MaskPlanes.insert(strIntPair{"SAT",   1}); /* Pixel is/was saturated */
    MaskPlanes.insert(strIntPair{"INTRP", 2}); /* Pixel has been interpolated */
    MaskPlanes.insert(strIntPair{"CR",    3}); /* Pixel is part of a cosmic ray */
    MaskPlanes.insert(strIntPair{"EDGE",  4}); /* Pixel is too close to the edge of the image, e.g. during convolution */
}
}
}
        
#endif // !defined(LSST_Ip_Diffim_MaskPlanes_H)
