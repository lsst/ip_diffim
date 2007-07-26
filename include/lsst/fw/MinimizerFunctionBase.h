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

#include <Minuit/FCNBase.h>
#include <lsst/fw/Function.h>

namespace lsst {
namespace fw {
namespace function {

    template<typename ReturnT>
    class MinimizerFunctionBase1 : public FCNBase {
        explicit MinimizerFunctionBase1();
        virtual ~MinimizerFunctionBase1() {};
        explicit MinimizerFunctionBase1(const std::vector<double>& measurementVector,
                                        const std::vector<double>& varianceVector,
                                        const std::vector<double>& positionVector,
                                        double errorDef,
                                        lsst::fw::function::Function1<ReturnT> theFunction
            );
        // Required by FCNBase
        virtual double up() const {return _errorDef;}
        virtual double operator() (const std::vector<double>&) const;
        
        inline std::vector<double> getMeasurements() const {return _measurementVector;}
        inline std::vector<double> getVariances() const {return _varianceVector;}
        inline std::vector<double> getPositions() const {return _positionVector;}
        inline void setErrorDef(double def) {_errorDef=def;}
    private:
        std::vector<double> _measurementVector;
        std::vector<double> _varianceVector;
        std::vector<double> _positionVector;
        lsst::fw::function::Function1<ReturnT> _theFunction;
        double _errorDef;
    };
        
    template<typename ReturnT>
    class MinimizerFunctionBase2 : public FCNBase {
        explicit MinimizerFunctionBase2();
        virtual ~MinimizerFunctionBase2() {};
        explicit MinimizerFunctionBase2(const std::vector<double>& measurementVector,
                                        const std::vector<double>& varianceVector,
                                        const std::vector<double>& position1Vector,
                                        const std::vector<double>& position2Vector,
                                        double errorDef,
                                        lsst::fw::function::Function2<ReturnT> theFunction
            );
        // Required by FCNBase
        virtual double up() const {return _errorDef;}
        virtual double operator() (const std::vector<double>&) const;
        
        inline std::vector<double> getMeasurements() const {return _measurementVector;}
        inline std::vector<double> getVariances() const {return _varianceVector;}
        inline std::vector<double> getPosition1() const {return _position1Vector;}
        inline std::vector<double> getPosition2() const {return _position2Vector;}
        inline void setErrorDef(double def) {_errorDef=def;}
    private:
        std::vector<double> _measurementVector;
        std::vector<double> _varianceVector;
        std::vector<double> _position1Vector;
        std::vector<double> _position2Vector;
        lsst::fw::function::Function2<ReturnT> _theFunction;
        double _errorDef;
    };
        
    
}   // namespace function        
}   // namespace fw
}   // namespace lsst

#include <lsst/fw/MinimizerFunctionBase.cc>

#endif // !defined(LSST_FW_Minimizer_H)
