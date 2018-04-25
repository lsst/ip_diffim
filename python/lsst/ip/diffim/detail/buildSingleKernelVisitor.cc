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
#include "pybind11/stl.h"

#include <memory>
#include <string>

#include <Eigen/Core>
#include "ndarray/pybind11.h"

#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
#include "lsst/pex/policy/Policy.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {
namespace detail {

namespace {

/**
 * Wrap BuildSingleKernelVisitor and the factory function for one pixel type
 *
 * @tparam PixelT  Pixel type, e.g. `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Class name suffix associated with PixeT, e.g. "F" for `float`
 */
template <typename PixelT>
void declareBuildSingleKernelVisitor(py::module& mod, std::string const& suffix) {
    py::class_<BuildSingleKernelVisitor<PixelT>, std::shared_ptr<BuildSingleKernelVisitor<PixelT>>,
               afw::math::CandidateVisitor>
            cls(mod, ("BuildSingleKernelVisitor" + suffix).c_str());

    cls.def(py::init<afw::math::KernelList, pex::policy::Policy>(), "basisList"_a, "policy"_a);
    cls.def(py::init<afw::math::KernelList, pex::policy::Policy, Eigen::MatrixXd const&>(), "basisList"_a,
            "policy"_a, "hMat"_a);

    cls.def("setSkipBuilt", &BuildSingleKernelVisitor<PixelT>::setSkipBuilt, "skip"_a);
    cls.def("getNRejected", &BuildSingleKernelVisitor<PixelT>::getNRejected);
    cls.def("getNProcessed", &BuildSingleKernelVisitor<PixelT>::getNProcessed);
    cls.def("reset", &BuildSingleKernelVisitor<PixelT>::reset);
    cls.def("processCandidate", &BuildSingleKernelVisitor<PixelT>::processCandidate, "candidate"_a);

    mod.def("makeBuildSingleKernelVisitor",
            (std::shared_ptr<BuildSingleKernelVisitor<PixelT>>(*)(afw::math::KernelList const&,
                                                                  pex::policy::Policy const&)) &
                    makeBuildSingleKernelVisitor<PixelT>,
            "basisList"_a, "policy"_a);
    mod.def("makeBuildSingleKernelVisitor",
            (std::shared_ptr<BuildSingleKernelVisitor<PixelT>>(*)(
                    afw::math::KernelList const&, pex::policy::Policy const&, Eigen::MatrixXd const&)) &
                    makeBuildSingleKernelVisitor<PixelT>,
            "basisList"_a, "policy"_a, "hMat"_a);
}

}  // namespace lsst::ip::diffim::detail::<anonymous>

PYBIND11_PLUGIN(buildSingleKernelVisitor) {
    py::module::import("lsst.afw.math");
    py::module::import("lsst.pex.policy");

    py::module mod("buildSingleKernelVisitor");

    declareBuildSingleKernelVisitor<float>(mod, "F");

    return mod.ptr();
}

}  // detail
}  // diffim
}  // ip
}  // lsst
