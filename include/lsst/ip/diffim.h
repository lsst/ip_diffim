/**
 * \file
 * \brief An include file to include the header files for lsst::ip::diffim
 */
#ifndef LSST_IP_DIFFIM_H
#define LSST_IP_DIFFIM_H

#include "lsst/ip/diffim/ImageSubtract.h"
#include "lsst/ip/diffim/BasisSets.h"
#include "lsst/ip/diffim/PsfMatchingFunctor.h"
#include "lsst/ip/diffim/SpatialModelKernel.h"
#include "lsst/ip/diffim/KernelCandidate.h"
#include "lsst/ip/diffim/KernelCandidateDetection.h"

#include "lsst/ip/diffim/KernelSumVisitor.h"
#include "lsst/ip/diffim/BuildSingleKernelVisitor.h"
#include "lsst/ip/diffim/KernelPcaVisitor.h"
#include "lsst/ip/diffim/BuildSpatialKernelVisitor.h"
#include "lsst/ip/diffim/AssessSpatialKernelVisitor.h"

#endif // LSST_IP_DIFFIM_H
