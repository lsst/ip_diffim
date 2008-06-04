// -*- LSST-C++ -*-
#ifndef LSST_FW_kriging_H
#define LSST_FW_kriging_H
/**
 * \file
 *
 * Implentation of spatial modelling
 *
 * \author Andrew Becker 
 *
 * \ingroup diffim
 */
#include <vario.h>  
#include <glvars.h>

namespace lsst{
namespace ip{
namespace diffim{

//    /* copied from glvars.h */
//    typedef enum {
//        NSP = 0,   /* initial value */
//        UIF,       /* variogram modelling user interface */
//        OKR, UKR, SKR, /* ordinary, universal or simple kriging */
//        IDW,       /* inverse distance interpolation */
//        MED,       /* (local) sample median or quantile */
//        NRS,       /* neighbourhood size */
//        LSLM,      /* uncorrelated (or weighted) linear model */
//        GSI, ISI,  /* Gaussian/indicator (conditional) simulation */
//        MAPVALUE,  /* mask map value at data location */
//        SEM, COV,  /* sample (cross) semivariance or covariance */
//        SPREAD,    /* distance to nearest sample */
//        XYP,       /* x and y coordinate of location */
//        POLY,      /* point-in-polygon */
//        DIV,       /* diversity, range */
//        SKEW,      /* skewness, kurtosis */
//        LSEM,      /* locally fitted semivariogram parameters */
//        TEST       /* does nothing really */
//    } METHOD;
    
    class Variogram : public lsst::mwi::data::LsstBase {
    public:
        explicit Variogram();
        virtual ~Variogram() {};

        /* gstat objects */
        VARIOGRAM *vp = NULL;
        DATA **d = NULL, *d1 = NULL, *d2 = NULL;

        fillVariogram();
        calcVariogram();
        doVariogram(int nvars, METHOD m);
    }
    
    class SpatialInterpolator : public lsst::mwi::data::LsstBase {
    public:
        explicit SpatialInteroplator();
        virtual ~SpatialInterpolator() {};
        
        void setMethod(METHOD);
        METHOD getMethod(void);

        Variogram *getVariogram(int i);

    private:
        
    };
} //namespace lsst
} //namespace ip
} //namespace diffim
        
