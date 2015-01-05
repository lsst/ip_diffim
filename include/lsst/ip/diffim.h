/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
/**
 * \file
 * \brief An include file to include the header files for lsst::ip::diffim
 */
#ifndef LSST_IP_DIFFIM_H
#define LSST_IP_DIFFIM_H

#include "lsst/ip/diffim/BasisLists.h"
#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/ImageStatistics.h"
#include "lsst/ip/diffim/FindSetBits.h"

#include "lsst/ip/diffim/KernelSolution.h"
#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelCandidateDetection.h"

#include "lsst/ip/diffim/KernelPca.h"
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
#include "lsst/ip/diffim/KernelSumVisitor.h"

#include "lsst/ip/diffim/DipoleAlgorithms.h"
#include "lsst/ip/diffim/SillyCentroid.h"

#endif // LSST_IP_DIFFIM_H
