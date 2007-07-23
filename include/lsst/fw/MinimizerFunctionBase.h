// -*- LSST-C++ -*-
#ifndef LSST_FW_Minimizer_H
#define LSST_FW_Minimizer_H
/**
 * \file
 *
 * Class that Minuit knows how to minimize, that contains an lsst::fw::function::Function
 *
 * \author Andrew Becker
 *
 * \ingroup fw
 */

#include "Minuit/FCNBase.h"
#include <lssw/fw/Function.h>

namespace lsst {
namespace fw {
namespace function {

    template<typename ReturnT>
    class MinimizerFunctionBase : public FCNBase {
        
    public:
        explicit MinimizerFunctionBase();
        virtual ~MinimizerFunctionBase() {};
        explicit MinimizerFunctionBase(const std::vector<double>& measurementVector,
                                       const std::vector<double>& positionVector,
                                       const std::vector<double>& varianceVector,
                                       double errorDef,
                                       lsst::fw::function::Function<ReturnT> theFunction
            );
        
        // Required by FCNBase
        virtual double up() const;
        virtual double operator() (const std::vector<double>&) const;
        
        virtual std::vector<double> getMeasurements() const {return _measurementVector;}
        virtual std::vector<double> getPositions() const {return _positionVector;}
        virtual std::vector<double> getVariances() const {return _varianceVector;}
        
        // For chi2 1-sigma, def = 1
        // For chi2 2-sigma, def = 4
        // Etc...
        inline void setErrorDef(double def);
        
    private:
        std::vector<double> _measurementVector;
        std::vector<double> _positionVector;
        std::vector<double> _varianceVector;
        double _errorDef;
        lsst::fw::function::Function<ReturnT> _theFunction;
    };
    
}   // namespace function        
}   // namespace fw
}   // namespace lsst

#include <lsst/fw/MinimizerFunctionBase.cc>

#endif // !defined(LSST_FW_Minimizer_H)
