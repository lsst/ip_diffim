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

#include <Eigen/Core>

#include "ndarray/pybind11.h"

#include "lsst/ip/diffim/BasisLists.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {

PYBIND11_PLUGIN(basisLists) {
    py::module::import("lsst.afw.math");
    py::module::import("lsst.pex.policy");

    py::module mod("basisLists");

    mod.def("makeDeltaFunctionBasisList", &makeDeltaFunctionBasisList, "width"_a, "height"_a);
    mod.def("makeRegularizationMatrix", &makeRegularizationMatrix, "policy"_a);
    mod.def("makeForwardDifferenceMatrix", &makeForwardDifferenceMatrix, "width"_a, "height"_a, "orders"_a,
            "borderPenalty"_a, "fitForBackground"_a);
    mod.def("makeCentralDifferenceMatrix", &makeCentralDifferenceMatrix, "width"_a, "height"_a, "stencil"_a,
            "borderPenalty"_a, "fitForBackground"_a);
    mod.def("renormalizeKernelList", &renormalizeKernelList, "kernelListIn"_a);
    mod.def("makeAlardLuptonBasisList", &makeAlardLuptonBasisList, "halfWidth"_a, "nGauss"_a, "sigGauss"_a,
            "degGauss"_a);

    return mod.ptr();
}

}  // diffim
}  // ip
}  // lsst
