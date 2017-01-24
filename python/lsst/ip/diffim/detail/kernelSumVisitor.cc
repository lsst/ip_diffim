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
#include "lsst/ip/diffim/KernelSumVisitor.h"
#include "lsst/pex/policy/Policy.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {
namespace detail {

namespace {

/**
 * Wrap KernelSumVisitor and the factory function for one pixel type
 *
 * @tparam PixelT  Pixel type, e.g. `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Class name suffix associated with PixeT, e.g. "F" for `float`
 */
template <typename PixelT>
void declareKernelSumVisitor(py::module& mod, std::string const& suffix) {

    using Class = KernelSumVisitor<PixelT>;

    py::class_<Class, std::shared_ptr<Class>, afw::math::CandidateVisitor>
        cls(mod, ("KernelSumVisitor" + suffix).c_str());

    py::enum_<typename Class::Mode>(cls, "Mode")
        .value("AGGREGATE", Class::Mode::AGGREGATE)
        .value("REJECT", Class::Mode::REJECT)
        .export_values();

    cls.def(py::init<pex::policy::Policy>(), "policy"_a);

    cls.def("setMode", &Class::setMode, "mode"_a);
    cls.def("getNRejected", &Class::getNRejected);
    cls.def("getkSumMean", &Class::getkSumMean);
    cls.def("getkSumStd", &Class::getkSumStd);
    cls.def("getdkSumMax", &Class::getdkSumMax);
    cls.def("getkSumNpts", &Class::getkSumNpts);
    cls.def("resetKernelSum", &Class::resetKernelSum);
    cls.def("processCandidate", &Class::processCandidate, "candidate"_a);
    cls.def("processKsumDistribution", &Class::processKsumDistribution);

    mod.def("makeKernelSumVisitor", &makeKernelSumVisitor<PixelT>, "policy"_a);
}

}  // namespace lsst::ip::diffim::detail::<anonymous>

PYBIND11_PLUGIN(_kernelSumVisitor) {
    py::module mod("_kernelSumVisitor", "Python wrapper for KernelSumVisitor.h");

    declareKernelSumVisitor<float>(mod, "F");

    return mod.ptr();
}
}
}
}
}  // lsst::ip::diffim::detail
