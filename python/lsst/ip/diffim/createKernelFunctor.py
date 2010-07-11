# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.pex.logging as pexLogging
import diffimLib 

def createKernelFunctor(policy):
    kCols = policy.getInt("kernelCols")
    kRows = policy.getInt("kernelRows")
    kType = policy.getString("kernelBasisSet")
    useRegularization = policy.getBool("useRegularization")
    
    if kType == "alard-lupton":
        alardNGauss   = policy.getInt("alardNGauss")
        alardSigGauss = policy.getDoubleArray("alardSigGauss")
        alardDegGauss = policy.getIntArray("alardDegGauss")

        if useRegularization:
            pexLogging.Trace("lsst.ip.diffim.createKernelFunctor", 1,
                             "Warning : Regularization not enabled for Alard-Lupton kernels")

        try:
            assert len(alardSigGauss) == alardNGauss
            assert len(alardDegGauss) == alardNGauss
            assert kCols == kRows          # square
            assert kCols % 2 == 1          # odd sized
        except:
            raise RuntimeError("Invalid Alard-Lupton kernel configuration")
        
        kHalfWidth = kCols//2
        basisList  = diffimLib.generateAlardLuptonBasisSet(kHalfWidth,
                                                           alardNGauss,
                                                           alardSigGauss,
                                                           alardDegGauss)
        kFunctor   = diffimLib.PsfMatchingFunctorF(basisList)
        return kFunctor

    elif kType == "delta-function":
        basisList = diffimLib.generateDeltaFunctionBasisSet(kCols, kRows)
        
        if useRegularization:
            regularizationOrder      = policy.getInt("regularizationOrder")
            regularizationBoundary   = policy.getInt("regularizationBoundary")
            regularizationDifference = policy.getInt("regularizationDifference")
            h = diffimLib.generateFiniteDifferenceRegularization(kCols, kRows,
                                                                 regularizationOrder,
                                                                 regularizationBoundary,
                                                                 regularizationDifference)

            kFunctor = diffimLib.PsfMatchingFunctorF(basisList, h)
            return kFunctor

        kFunctor = diffimLib.PsfMatchingFunctorF(basisList)
        return kFunctor

    else:
        raise RuntimeError("Invalid kernel type : %s" % (kType))
    
