/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"

#include <memory>

#include "Eigen/Core"
#include "ndarray/pybind11.h"

#include "lsst/pex/policy/Policy.h"
#include "lsst/ip/diffim/KernelSolution.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {

namespace {

/**
 * Wrap KernelSolution
 */
void declareKernelSolution(py::module &mod) {
    py::class_<KernelSolution, std::shared_ptr<KernelSolution>> cls(mod, "KernelSolution");

    cls.def(py::init<Eigen::MatrixXd, Eigen::VectorXd, bool>(), "mMat"_a, "bVec"_a, "fitForBackground"_a);
    cls.def(py::init<bool>(), "fitForBackground"_a);
    cls.def(py::init<>());

    py::enum_<KernelSolution::KernelSolvedBy>(cls, "KernelSolvedBy")
            .value("NONE", KernelSolution::KernelSolvedBy::NONE)
            .value("CHOLESKY_LDLT", KernelSolution::KernelSolvedBy::CHOLESKY_LDLT)
            .value("CHOLESKY_LLT", KernelSolution::KernelSolvedBy::CHOLESKY_LLT)
            .value("LU", KernelSolution::KernelSolvedBy::LU)
            .value("EIGENVECTOR", KernelSolution::KernelSolvedBy::EIGENVECTOR)
            .export_values();

    py::enum_<KernelSolution::ConditionNumberType>(cls, "ConditionNumberType")
            .value("EIGENVALUE", KernelSolution::ConditionNumberType::EIGENVALUE)
            .value("SVD", KernelSolution::ConditionNumberType::SVD)
            .export_values();

    cls.def("solve", (void (KernelSolution::*)()) & KernelSolution::solve);
    cls.def("solve", (void (KernelSolution::*)(Eigen::MatrixXd const &, Eigen::VectorXd const &)) &
                             KernelSolution::solve,
            "mMat"_a, "bVec"_a);
    cls.def("getSolvedBy", &KernelSolution::getSolvedBy);
    cls.def("getConditionNumber", (double (KernelSolution::*)(KernelSolution::ConditionNumberType)) &
                                          KernelSolution::getConditionNumber,
            "conditionType"_a);
    cls.def("getConditionNumber",
            (double (KernelSolution::*)(Eigen::MatrixXd const &, KernelSolution::ConditionNumberType)) &
                    KernelSolution::getConditionNumber,
            "mMat"_a, "conditionType"_a);
    cls.def("getM", &KernelSolution::getM, py::return_value_policy::copy);
    cls.def("getB", &KernelSolution::getB, py::return_value_policy::copy);
    cls.def("printM", &KernelSolution::printM);
    cls.def("printB", &KernelSolution::printB);
    cls.def("printA", &KernelSolution::printA);
    cls.def("getId", &KernelSolution::getId);
}

/**
 * Wrap StaticKernelSolution
 *
 * @tparam InputT  Pixel type of science and template images, e,g, `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Suffix associated with `InputT`, e.g. "F" for `float`
 */
template <typename InputT>
void declareStaticKernelSolution(py::module &mod, std::string const &suffix) {
    py::class_<StaticKernelSolution<InputT>, std::shared_ptr<StaticKernelSolution<InputT>>, KernelSolution>
            cls(mod, ("StaticKernelSolution" + suffix).c_str());

    cls.def(py::init<lsst::afw::math::KernelList const &, bool>(), "basisList"_a, "fitForBackground"_a);

    cls.def("solve", (void (StaticKernelSolution<InputT>::*)()) & StaticKernelSolution<InputT>::solve);
    cls.def("build", &StaticKernelSolution<InputT>::build, "templateImage"_a, "scienceImage"_a,
            "varianceEstimate"_a);
    cls.def("getKernel", &StaticKernelSolution<InputT>::getKernel);
    cls.def("makeKernelImage", &StaticKernelSolution<InputT>::makeKernelImage);
    cls.def("getBackground", &StaticKernelSolution<InputT>::getBackground);
    cls.def("getKsum", &StaticKernelSolution<InputT>::getKsum);
    cls.def("getSolutionPair", &StaticKernelSolution<InputT>::getSolutionPair);
}

/**
 * Wrap MaskedKernelSolution
 *
 * @tparam InputT  Pixel type of science and template images, e,g, `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Suffix associated with `InputT`, e.g. "F" for `float`
 */
template <typename InputT>
void declareMaskedKernelSolution(py::module &mod, std::string const &suffix) {
    py::class_<MaskedKernelSolution<InputT>, std::shared_ptr<MaskedKernelSolution<InputT>>,
               StaticKernelSolution<InputT>>
            cls(mod, ("MaskedKernelSolution" + suffix).c_str());

    cls.def(py::init<lsst::afw::math::KernelList const &, bool>(), "basisList"_a, "fitForBackground"_a);

    cls.def("buildOrig", &MaskedKernelSolution<InputT>::buildOrig, "templateImage"_a, "scienceImage"_a,
            "varianceEstimate"_a, "pixelMask"_a);
    cls.def("buildWithMask", &MaskedKernelSolution<InputT>::buildWithMask, "templateImage"_a,
            "scienceImage"_a, "varianceEstimate"_a, "pixelMask"_a);
    cls.def("buildSingleMaskOrig", &MaskedKernelSolution<InputT>::buildSingleMaskOrig, "templateImage"_a,
            "scienceImage"_a, "varianceEstimate"_a, "maskBox"_a);
}

/**
 * Wrap RegularizedKernelSolution
 *
 * @tparam InputT  Pixel type of science and template images, e,g, `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Suffix associated with `InputT`, e.g. "F" for `float`
 */
template <typename InputT>
void declareRegularizedKernelSolution(py::module &mod, std::string const &suffix) {
    py::class_<RegularizedKernelSolution<InputT>, std::shared_ptr<RegularizedKernelSolution<InputT>>,
               StaticKernelSolution<InputT>>
            cls(mod, ("RegularizedKernelSolution" + suffix).c_str());

    cls.def(py::init<lsst::afw::math::KernelList const &, bool, Eigen::MatrixXd const &,
                     pex::policy::Policy>(),
            "basisList"_a, "fitForBackground"_a, "hMat"_a, "policy"_a);

    cls.def("solve",
            (void (RegularizedKernelSolution<InputT>::*)()) & RegularizedKernelSolution<InputT>::solve);
    cls.def("getLambda", &RegularizedKernelSolution<InputT>::getLambda);
    cls.def("estimateRisk", &RegularizedKernelSolution<InputT>::estimateRisk, "maxCond"_a);
    cls.def("getM", &RegularizedKernelSolution<InputT>::getM);
}

/**
 * Wrap SpatialKernelSolution
 */
void declareSpatialKernelSolution(py::module &mod) {
    py::class_<SpatialKernelSolution, std::shared_ptr<SpatialKernelSolution>, KernelSolution> cls(
            mod, "SpatialKernelSolution");

    cls.def(py::init<lsst::afw::math::KernelList const &, lsst::afw::math::Kernel::SpatialFunctionPtr,
                     lsst::afw::math::Kernel::SpatialFunctionPtr, pex::policy::Policy>(),
            "basisList"_a, "spatialKernelFunction"_a, "background"_a, "policy"_a);

    cls.def("solve", (void (SpatialKernelSolution::*)()) & SpatialKernelSolution::solve);
    cls.def("addConstraint", &SpatialKernelSolution::addConstraint, "xCenter"_a, "yCenter"_a, "qMat"_a,
            "wVec"_a);
    cls.def("makeKernelImage", &SpatialKernelSolution::makeKernelImage, "pos"_a);
    cls.def("getSolutionPair", &SpatialKernelSolution::getSolutionPair);
}

}  // namespace lsst::ip::diffim::<anonymous>

PYBIND11_PLUGIN(kernelSolution) {
    py::module::import("lsst.afw.geom");
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.math");
    py::module::import("lsst.pex.policy");

    py::module mod("kernelSolution");

    declareKernelSolution(mod);
    declareStaticKernelSolution<float>(mod, "F");
    declareMaskedKernelSolution<float>(mod, "F");
    declareRegularizedKernelSolution<float>(mod, "F");
    declareSpatialKernelSolution(mod);

    return mod.ptr();
}

}  // diffim
}  // ip
}  // lsst
