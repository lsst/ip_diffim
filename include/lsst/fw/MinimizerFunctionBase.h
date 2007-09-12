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
#include <lsst/mwi/data/LsstData.h>
#include <lsst/fw/Function.h>
#include <boost/shared_ptr.hpp>

namespace lsst {
namespace fw {
namespace function {

    template<typename ReturnT>
    class MinimizerFunctionBase1 : public FCNBase, private lsst::mwi::data::LsstBase {
    public:
        explicit MinimizerFunctionBase1();
        virtual ~MinimizerFunctionBase1() {};
        explicit MinimizerFunctionBase1(std::vector<double> const &measurementVector,
                                        std::vector<double> const &varianceVector,
                                        std::vector<double> const &positionVector,
                                        double errorDef,
                                        boost::shared_ptr<lsst::fw::function::Function1<ReturnT> > theFunctionPtr
            );
        // Required by FCNBase
        virtual double up() const {return _errorDef;}
        virtual double operator() (const std::vector<double>&) const;

        //void minimizee(std::vector<double> &parameters,
        //std::vector<std::pair<double,double> > &errors);
        inline std::vector<double> getMeasurements() const {return _measurementVector;}
        inline std::vector<double> getVariances() const {return _varianceVector;}
        inline std::vector<double> getPositions() const {return _positionVector;}
        inline void setErrorDef(double def) {_errorDef=def;}
    private:
        std::vector<double> _measurementVector;
        std::vector<double> _varianceVector;
        std::vector<double> _positionVector;
        double _errorDef;
        boost::shared_ptr<lsst::fw::function::Function1<ReturnT> > _theFunctionPtr;
    };
        
    template<typename ReturnT>
    class MinimizerFunctionBase2 : public FCNBase, private lsst::mwi::data::LsstBase {
    public:
        explicit MinimizerFunctionBase2();
        virtual ~MinimizerFunctionBase2() {};
        explicit MinimizerFunctionBase2(std::vector<double> const &measurementVector,
                                        std::vector<double> const &varianceVector,
                                        std::vector<double> const &position1Vector,
                                        std::vector<double> const &position2Vector,
                                        double errorDef,
                                        boost::shared_ptr<lsst::fw::function::Function2<ReturnT> > theFunctionPtr
            );
        // Required by FCNBase
        virtual double up() const {return _errorDef;}
        virtual double operator() (const std::vector<double>&) const;
        
        //void minimizee(std::vector<double> &parameters,
        //std::vector<std::pair<double,double> > &errors);
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
        double _errorDef;
        boost::shared_ptr<lsst::fw::function::Function2<ReturnT> > _theFunctionPtr;
    };
        
    template<typename ReturnT>
    void minimize(
        lsst::fw::function::MinimizerFunctionBase1<ReturnT> &theFunction,
        std::vector<double> &parameters,
        std::vector<double> const &stepsize,
        std::vector<std::pair<double,double> > &errors
        );

    template<typename ReturnT>
    void minimize(
        lsst::fw::function::MinimizerFunctionBase2<ReturnT> &theFunction,
        std::vector<double> &parameters,
        std::vector<double> const &stepsize,
        std::vector<std::pair<double,double> > &errors
        );
        
    
}   // namespace function        
}   // namespace fw
}   // namespace lsst

#include <lsst/fw/MinimizerFunctionBase.cc>

#endif // !defined(LSST_FW_Minimizer_H)
