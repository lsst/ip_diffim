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

#include "lsst/ip/diffim/DipoleAlgorithms.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/pex/config/python.h"  // for LSST_DECLARE_CONTROL_FIELD

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace ip {
namespace diffim {

namespace {

void declareDipoleCentroidControl(py::module &mod) {
    py::class_<DipoleCentroidControl, std::shared_ptr<DipoleCentroidControl>> cls(mod,
                                                                                  "DipoleCentroidControl");

    cls.def(py::init<>());
}

void declareDipoleFluxControl(py::module &mod) {
    py::class_<DipoleFluxControl, std::shared_ptr<DipoleFluxControl>> cls(mod, "DipoleFluxControl");

    cls.def(py::init<>());
}

void declareDipolePsfFluxControl(py::module &mod) {
    py::class_<PsfDipoleFluxControl, std::shared_ptr<PsfDipoleFluxControl>, DipoleFluxControl> cls(
            mod, "PsfDipoleFluxControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, PsfDipoleFluxControl, stepSizeCoord);
    LSST_DECLARE_CONTROL_FIELD(cls, PsfDipoleFluxControl, stepSizeFlux);
    LSST_DECLARE_CONTROL_FIELD(cls, PsfDipoleFluxControl, errorDef);
    LSST_DECLARE_CONTROL_FIELD(cls, PsfDipoleFluxControl, maxFnCalls);
}

void declareDipoleCentroidAlgorithm(py::module &mod) {
    // Abstract class, so add a leading underscore to Python name and do not wrap constructor
    py::class_<DipoleCentroidAlgorithm, std::shared_ptr<DipoleCentroidAlgorithm>, meas::base::SimpleAlgorithm>
            cls(mod, "_DipoleCentroidAlgorithm");

    cls.attr("FAILURE") = py::cast(DipoleCentroidAlgorithm::FAILURE);
    cls.attr("POS_FLAG") = py::cast(DipoleCentroidAlgorithm::POS_FLAG);
    cls.attr("NEG_FLAG") = py::cast(DipoleCentroidAlgorithm::NEG_FLAG);
    cls.def_static("getFlagDefinitions", &DipoleCentroidAlgorithm::getFlagDefinitions,
                   py::return_value_policy::copy);

    cls.def("getPositiveKeys", &DipoleCentroidAlgorithm::getPositiveKeys);
    cls.def("getNegativeKeys", &DipoleCentroidAlgorithm::getNegativeKeys);
}

void declareDipoleFluxAlgorithm(py::module &mod) {
    // Abstract class, so add a leading underscore to Python name and do not wrap constructor
    py::class_<DipoleFluxAlgorithm, std::shared_ptr<DipoleFluxAlgorithm>, meas::base::SimpleAlgorithm> cls(
            mod, "_DipoleFluxAlgorithm");

    cls.attr("FAILURE") = py::cast(DipoleFluxAlgorithm::FAILURE);
    cls.attr("POS_FLAG") = py::cast(DipoleFluxAlgorithm::POS_FLAG);
    cls.attr("NEG_FLAG") = py::cast(DipoleFluxAlgorithm::NEG_FLAG);
    cls.def_static("getFlagDefinitions", &DipoleFluxAlgorithm::getFlagDefinitions,
                   py::return_value_policy::copy);

    cls.def("getPositiveKeys", &DipoleFluxAlgorithm::getPositiveKeys);
    cls.def("getNegativeKeys", &DipoleFluxAlgorithm::getNegativeKeys);
}

void declareNaiveDipoleFlux(py::module &mod) {
    py::class_<NaiveDipoleFlux, std::shared_ptr<NaiveDipoleFlux>, DipoleFluxAlgorithm> cls(mod,
                                                                                           "NaiveDipoleFlux");

    cls.def(py::init<NaiveDipoleFlux::Control const &, std::string const &, afw::table::Schema &>(), "ctrl"_a,
            "name"_a, "schema"_a);

    cls.def("measure", &NaiveDipoleFlux::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &NaiveDipoleFlux::fail, "measRecord"_a, "error"_a = NULL);
}

void declareNaiveDipoleCentroid(py::module &mod) {
    py::class_<NaiveDipoleCentroid, std::shared_ptr<NaiveDipoleCentroid>, DipoleCentroidAlgorithm> cls(
            mod, "NaiveDipoleCentroid");

    cls.def(py::init<NaiveDipoleCentroid::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);

    cls.def("getCenterKeys", &NaiveDipoleCentroid::getCenterKeys);
    cls.def("getPositiveKeys", &NaiveDipoleCentroid::getPositiveKeys);
    cls.def("getNegativeKeys", &NaiveDipoleCentroid::getNegativeKeys);

    cls.def("measure", &NaiveDipoleCentroid::measure, "measRecord"_a, "exposure"_a);
    cls.def("mergeCentroids", &NaiveDipoleCentroid::mergeCentroids, "source"_a, "posValue"_a, "negValue"_a);
    cls.def("fail", &NaiveDipoleCentroid::fail, "measRecord"_a, "error"_a = NULL);
}

void declarePsfDipoleFlux(py::module &mod) {
    py::class_<PsfDipoleFlux, std::shared_ptr<PsfDipoleFlux>, DipoleFluxAlgorithm> cls(mod, "PsfDipoleFlux");

    cls.def(py::init<PsfDipoleFlux::Control const &, std::string const &, afw::table::Schema &>(), "ctrl"_a,
            "name"_a, "schema"_a);

    cls.def("chi2", &PsfDipoleFlux::chi2, "source"_a, "exposure"_a, "negCenterX"_a, "negCenterY"_a,
            "negFlux"_a, "posCenterX"_a, "posCenterY"_a, "posFlux"_a);
    cls.def("measure", &PsfDipoleFlux::measure, "measRecord"_a, "exposure"_a);
    cls.def("fail", &PsfDipoleFlux::fail, "measRecord"_a, "error"_a = NULL);
}

}  // namespace lsst::ip::diffim::<anonymous>

PYBIND11_PLUGIN(_dipoleAlgorithms) {
    py::module mod("_dipoleAlgorithms", "Python wrapper for DipoleAlgorithms.h");

    declareDipoleCentroidControl(mod);
    declareDipoleFluxControl(mod);
    declareDipolePsfFluxControl(mod);
    declareDipoleCentroidAlgorithm(mod);
    declareDipoleFluxAlgorithm(mod);
    declareNaiveDipoleFlux(mod);
    declareNaiveDipoleCentroid(mod);
    declarePsfDipoleFlux(mod);

    return mod.ptr();
}

}  // diffim
}  // ip
}  // lsst
