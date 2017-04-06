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

#include "lsst/ip/diffim/KernelCandidateDetection.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {

namespace {

/**
 * Wrap `KernelCandidateDetection` for one pixel type
 *
 * @tparam PixelT  Pixel type of image plane of MaskedImage, e.g. `float`
 * @param mod  pybind11 module
 * @param[in] suffix  Class name suffix associated with PixeT, e.g. "F" for `float`
 */
template <typename PixelT>
void declareKernelCandidateDetection(py::module &mod, std::string const &suffix) {
    py::class_<KernelCandidateDetection<PixelT>, std::shared_ptr<KernelCandidateDetection<PixelT>>> cls(
            mod, ("KernelCandidateDetection" + suffix).c_str());

    cls.def(py::init<pex::policy::Policy const &>(), "policy"_a);

    cls.def("apply", &KernelCandidateDetection<PixelT>::apply, "templateMaskedImage"_a,
            "scienceMaskedImage"_a);
    cls.def("growCandidate", &KernelCandidateDetection<PixelT>::growCandidate, "footprint"_a, "fpGrowPix"_a,
            "templateMaskedImage"_a, "scienceMaskedImage"_a);
    cls.def("getFootprints", &KernelCandidateDetection<PixelT>::getFootprints);
}

}  // namespace lsst::ip::diffim::<anonymous>

PYBIND11_PLUGIN(kernelCandidateDetection) {
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.detection");
    py::module::import("lsst.pex.policy");

    py::module mod("kernelCandidateDetection");

    declareKernelCandidateDetection<float>(mod, "F");

    return mod.ptr();
}

}  // diffim
}  // ip
}  // lsst
