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
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>

#include "lsst/afw/math/SpatialCell.h"
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {
namespace detail {

namespace {

/**
 * Wrap AssessSpatialKernelVisitor and the factory function for one pixel type
 *
 * @tparam PixelT  Pixel type, e.g. `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Class name suffix associated with PixeT, e.g. "F" for `float`
 */
template <typename PixelT>
void declareAssessSpatialKernelVisitor(py::module& mod, std::string const& suffix) {
    py::class_<AssessSpatialKernelVisitor<PixelT>, std::shared_ptr<AssessSpatialKernelVisitor<PixelT>>,
               lsst::afw::math::CandidateVisitor>
            cls(mod, ("AssessSpatialKernelVisitor" + suffix).c_str());

    cls.def(py::init<std::shared_ptr<afw::math::LinearCombinationKernel>,
                     afw::math::Kernel::SpatialFunctionPtr, pex::policy::Policy const&>(),
            "spatialKernel"_a, "spatialBackground"_a, "policy"_a);

    cls.def("reset", &AssessSpatialKernelVisitor<PixelT>::reset);
    cls.def("getNGood", &AssessSpatialKernelVisitor<PixelT>::getNGood);
    cls.def("getNRejected", &AssessSpatialKernelVisitor<PixelT>::getNRejected);
    cls.def("getNProcessed", &AssessSpatialKernelVisitor<PixelT>::getNProcessed);
    cls.def("processCandidate", &AssessSpatialKernelVisitor<PixelT>::processCandidate, "candidate"_a);

    mod.def("makeAssessSpatialKernelVisitor", &makeAssessSpatialKernelVisitor<PixelT>, "spatialKernel"_a,
            "spatialBackground"_a, "policy"_a);
}

}  // namespace lsst::ip::diffim::detail::<anonymous>

PYBIND11_PLUGIN(_assessSpatialKernelVisitor) {
    py::module mod("_assessSpatialKernelVisitor", "Python wrapper for AssessSpatialKernelVisitor.h");

    declareAssessSpatialKernelVisitor<float>(mod, "F");

    return mod.ptr();
}
}
}
}
}  // lsst::ip::diffim::detail
