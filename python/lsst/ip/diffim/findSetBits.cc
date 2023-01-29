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
#include "lsst/cpputils/python.h"

#include <string>

#include "lsst/afw/image/LsstImageTypes.h"
#include "lsst/ip/diffim/FindSetBits.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {

namespace {

/**
 * Wrap FindSetBits for one mask pixel type
 *
 * @tparam MaskT  Mask type, typically lsst::afw::image::Mask<lsst::afw::image::MaskPixel>
 * @param mod  pybind11 module
 * @param[in] suffix  Class name suffix associated with mask pixel type, use "U" for `afw::image::MaskPixel`
 */
template <typename MaskT>
void declareFindSetBits(lsst::cpputils::python::WrapperCollection &wrappers, std::string const& suffix) {
    using PyFindSetBits = py::class_<FindSetBits<MaskT>>;

    std::string name = "FindSetBits" + suffix;
    wrappers.wrapType(PyFindSetBits(wrappers.module, name.c_str()), [](auto &mod, auto &cls) {

        cls.def(py::init<>());

        cls.def("reset", &FindSetBits<MaskT>::reset);
        cls.def("getBits", &FindSetBits<MaskT>::getBits);
        cls.def("apply", &FindSetBits<MaskT>::apply, "mask"_a);
    });
}

}  // namespace lsst::ip::diffim::<anonymous>

void wrapFindSetBits(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareFindSetBits<afw::image::Mask<afw::image::MaskPixel>>(wrappers, "U");
}

}  // diffim
}  // ip
}  // lsst
