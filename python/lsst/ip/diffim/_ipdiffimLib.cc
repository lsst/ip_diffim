/*
 * This file is part of ip_diffim.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"

namespace py = pybind11;
using namespace pybind11::literals;
using lsst::cpputils::python::WrapperCollection;

namespace lsst {
namespace ip {
namespace diffim {
void wrapBasisLists(WrapperCollection &wrappers);
void wrapDipoleAlgorithms(WrapperCollection &wrappers);
void wrapImageStatistics(WrapperCollection &wrappers);
void wrapImageSubtract(WrapperCollection &wrappers);
void wrapKernelCandidate(WrapperCollection &wrappers);
void wrapKernelSolution(WrapperCollection &wrappers);
namespace detail {
    void wrapAssessSpatialKernelVisitor(WrapperCollection &wrappers);
    void wrapBuildSingleKernelVisitor(WrapperCollection &wrappers);
    void wrapBuildSpatialKernelVisitor(WrapperCollection &wrappers);
    void wrapKernelPca(WrapperCollection &wrappers);
    void wrapKernelSumVisitor(WrapperCollection &wrappers);
}
PYBIND11_MODULE(_ipdiffimLib, mod) {
    lsst::utils::python::WrapperCollection wrappers(mod, "lsst.ip.diffim");
    wrappers.addInheritanceDependency("lsst.meas.base");
    wrappers.addInheritanceDependency("lsst.afw.math");
    wrappers.addSignatureDependency("lsst.afw.table");
    wrappers.addSignatureDependency("lsst.afw.image");
    wrappers.addSignatureDependency("lsst.afw.geom");
    wrappers.addSignatureDependency("lsst.daf.base");
    wrappers.addSignatureDependency("lsst.afw.detection");

    wrapBasisLists(wrappers);
    wrapDipoleAlgorithms(wrappers);
    wrapImageStatistics(wrappers);
    wrapImageSubtract(wrappers);
    wrapKernelCandidate(wrappers);
    wrapKernelSolution(wrappers);
    detail::wrapAssessSpatialKernelVisitor(wrappers);
    detail::wrapBuildSingleKernelVisitor(wrappers);
    detail::wrapBuildSpatialKernelVisitor(wrappers);
    detail::wrapKernelPca(wrappers);
    detail::wrapKernelSumVisitor(wrappers);
    wrappers.finish();
}
}
}
}
