// -*- LSST-C++ -*-
/**
 * \file
 *
 * Definition of member functions declared in MinimizerFunctionBase.h
 *
 * This file is meant to be included by lsst/fw/MinimizerFunctionBase.h
 *
 * \author Andrew Becker
 *
 * \ingroup fw
 */


// Constructors
template<typename ReturnT>
lsst::fw::function::MinimizerFunctionBase1<ReturnT>::MinimizerFunctionBase1()
:
    lsst::fw::LsstBase(typeid(this)),
    _measurementVector(),
    _varianceVector(),
    _positionVector(),
    _errorDef(1.),
    _theFunction()
{}

template<typename ReturnT>
lsst::fw::function::MinimizerFunctionBase1<ReturnT>::MinimizerFunctionBase1(
    const std::vector<double>& measurementVector,
    const std::vector<double>& varianceVector,
    const std::vector<double>& positionVector, 
    double errorDef,
    lsst::fw::function::Function1<ReturnT> theFunction)
:
    lsst::fw::LsstBase(typeid(this)),
    _measurementVector(measurementVector),
    _varianceVector(varianceVector),
    _positionVector(positionVector),
    _errorDef(errorDef),
    _theFunction(theFunction)
{}

template<typename ReturnT>
lsst::fw::function::MinimizerFunctionBase2<ReturnT>::MinimizerFunctionBase2()
:
    lsst::fw::LsstBase(typeid(this)),
    _measurementVector(),
    _varianceVector(),
    _position1Vector(),
    _position2Vector(),
    _errorDef(1.),
    _theFunction()
{}

template<typename ReturnT>
lsst::fw::function::MinimizerFunctionBase2<ReturnT>::MinimizerFunctionBase2(
    const std::vector<double>& measurementVector,
    const std::vector<double>& varianceVector,
    const std::vector<double>& position1Vector,
    const std::vector<double>& position2Vector,
    double errorDef,
    lsst::fw::function::Function2<ReturnT> theFunction)
:
    lsst::fw::LsstBase(typeid(this)),
    _measurementVector(measurementVector),
    _varianceVector(varianceVector),
    _position1Vector(position1Vector),
    _position2Vector(position2Vector),
    _errorDef(errorDef),
    _theFunction(theFunction)
{}



// Only method we need to set up; basically this is a chi^2 routine
template<typename ReturnT>
double lsst::fw::function::MinimizerFunctionBase1<ReturnT>::operator() (const std::vector<double>& par) const {
    // Initialize the function with the fit parameters
    this->_theFunction.setParameters(par);
    
    double chi2 = 0.;
    double resid;
    for (unsigned int i = 0; i < this->_theMeasurements.size(); i++) {
        resid = this->_theFunction(this->_positionVector[i]) - this->_measurementVector[i];
        chi2 += resid * resid / this->_varianceVector[i];
    }
    
    return chi2;
}


template<typename ReturnT>
double lsst::fw::function::MinimizerFunctionBase2<ReturnT>::operator() (const std::vector<double>& par) const {
    // Initialize the function with the fit parameters
    this->_theFunction.setParameters(par);
    
    double chi2 = 0.;
    double resid;
    for (unsigned int i = 0; i < this->_theMeasurements.size(); i++) {
        resid = this->_theFunction(this->_position1Vector[i], this->_position1Vector[i]) - this->_measurementVector[i];
        chi2 += resid * resid / this->_varianceVector[i];
    }
    
    return chi2;
}

