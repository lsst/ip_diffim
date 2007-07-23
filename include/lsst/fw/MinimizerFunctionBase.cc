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
lsst::fw::function::MinimizerFunctionBase()
:
    lsst::fw::function::LsstBase(typeid(this)),
    _measurementVector(),
    _positionVector(),
    _varianceVector(),
    _errorDef(1.),
    _theFunction()
{}

template<typename ReturnT>
lsst::fw::function::MinimizerFunctionBase(
    const std::vector<double>& measurementVector,
    const std::vector<double>& positionVector,
    const std::vector<double>& varianceVector,
    double errorDef,
    lsst::fw::function::Function<ReturnT> theFunction)
:
    lsst::fw::function::LsstBase(typeid(this)),
    _measurementVector(measurementVector),
    _positionVector(positionVector),
    _varianceVector(varianceVector),
    _errorDef(errorDef),
    _theFunction(theFunction)
{}


// Required by FCNBase
template<typename ReturnT>
inline double lsst::fw::function::MinimizerFunctionBase<ReturnT>::up() const {
    return _errorDef;
}

template<typename ReturnT>
inline double lsst::fw::function::MinimizerFunctionBase<ReturnT>::operator(const std::vector<double>&) const {
    
}        

// Others
template<typename ReturnT>
inline std::vector<double> lsst::fw::function::MinimizerFunctionBase<ReturnT>::getMeasurements() const {
    return _measurementVector;
}

template<typename ReturnT>
inline std::vector<double> lsst::fw::function::MinimizerFunctionBase<ReturnT>::getPositions() const {
    return _positionVector;
}

template<typename ReturnT>
inline std::vector<double> lsst::fw::function::MinimizerFunctionBase<ReturnT>::getVariances() const {
    return _varianceVector;
}

template<typename ReturnT>
inline void lsst::fw::function::MinimizerFunctionBase<ReturnT>::setErrorDef(double def) {
    this->_errorDef(def);
}
