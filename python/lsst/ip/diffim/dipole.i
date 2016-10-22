// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
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
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%{
#include "lsst/afw/geom.h"
#include "lsst/afw/math.h"
#include "lsst/afw/table.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/meas/base.h"
%}

%declareNumPyConverters(lsst::meas::base::CentroidCov);

%include "std_vector.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/table/tableLib.i"
%import "lsst/pex/config.h"

%{
#include "lsst/ip/diffim/DipoleAlgorithms.h"
%}

/******************************************************************************/

%feature("abstract") lsst::ip::diffim::DipoleCentroidAlgorithm;
%feature("abstract") lsst::ip::diffim::DipoleFluxAlgorithm;
%feature("notabstract") lsst::ip::diffim::PsfDipoleFlux;
%feature("notabstract") lsst::ip::diffim::NaiveDipoleFlux;
%feature("notabstract") lsst::ip::diffim::NaiveDipoleCentroid;


%shared_ptr(lsst::ip::diffim::DipoleCentroidControl)
%shared_ptr(lsst::ip::diffim::DipoleFluxControl)
%shared_ptr(lsst::ip::diffim::DipoleCentroidAlgorithm)
%shared_ptr(lsst::ip::diffim::DipoleCentroidControl)
%shared_ptr(lsst::ip::diffim::DipoleFluxAlgorithm)
%shared_ptr(lsst::ip::diffim::DipoleFluxControl)
%shared_ptr(lsst::ip::diffim::NaiveDipoleCentroidControl)
%shared_ptr(lsst::ip::diffim::NaiveDipoleFluxControl)
%shared_ptr(lsst::ip::diffim::PsfDipoleFluxControl)
%shared_ptr(lsst::ip::diffim::NaiveDipoleCentroid)
%shared_ptr(lsst::ip::diffim::PsfDipoleFlux)
%shared_ptr(lsst::ip::diffim::NaiveDipoleFlux)

%include "lsst/ip/diffim/DipoleAlgorithms.h"
