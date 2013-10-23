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

"""Quality Assessment class for Kernel Candidates"""
import numpy as np
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
from . import diffimLib

class KernelCandidateQa(object):

    def __init__(self, nKernelSpatial):
        """Class to undertake QA of KernelCandidates after modeling of
        the Psf-matching kernel.  Both directly--fitted diffim (LOCAL)
        and spatially--interpolated kernel diffim (SPATIAL) metrics
        are calculated, based on the distribution of residuals in the
        KernelCandidates stamp.

        @param nKernelSpatial : Number of terms in the spatial model; needed to initialize per-basis QA arrays
        """
        self.fields = []
        self.fields.append(afwTable.Field["PointD"]("RegisterRefPosition",
                                              "Position of reference object for registration (radians)."))
        #TODO check units of the following angles
        self.fields.append(afwTable.Field["Angle"]("RegisterResidualBearing",
                                              "Angle of residual wrt declination parallel in radians"))

        self.fields.append(afwTable.Field["Angle"]("RegisterResidualDistance",
                                              "Offset of residual in radians"))
        for kType in ("LOCAL", "SPATIAL"):
            self.fields.append(afwTable.Field["F"]("KCDiffimMean_%s"%(kType),
                                                   "Mean of KernelCandidate diffim", "sigma"))

            self.fields.append(afwTable.Field["F"]("KCDiffimMedian_%s"%(kType),
                                                   "Median of KernelCandidate diffim", "sigma"))

            self.fields.append(afwTable.Field["F"]("KCDiffimIQR_%s"%(kType),
                                                   "Inner quartile range of KernelCandidate diffim",
                                                   "sigma"))

            self.fields.append(afwTable.Field["F"]("KCDiffimStDev_%s"%(kType),
                                                   "Standard deviation of KernelCandidate diffim",
                                                   "sigma"))

            self.fields.append(afwTable.Field["F"]("KCDiffimKSD_%s"%(kType),
                                                   "D from K-S test of diffim pixels relative to Normal"))

            self.fields.append(afwTable.Field["F"]("KCDiffimKSProb_%s"%(kType),
                                                   "Prob from K-S test of diffim pixels relative to Normal",
                                                   "likelihood"))

            self.fields.append(afwTable.Field["F"]("KCDiffimADA2_%s"%(kType),
                    "Anderson-Darling test statistic of diffim pixels relative to Normal"))

            self.fields.append(afwTable.Field["ArrayD"]("KCDiffimADCrit_%s"%(kType),
                    "Critical values for the significance levels in KCDiffimADSig.  If A2 is greater than this number, hypothesis that the two distributions are the same can be rejected.", 5))

            self.fields.append(afwTable.Field["ArrayD"]("KCDiffimADSig_%s"%(kType),
                    "Anderson-Darling significance levels for the Normal distribution", 5))

            self.fields.append(afwTable.Field["F"]("KCDiffimChiSq_%s"%(kType),
                                                   "Reduced chi^2 of the residual.", "likelihood"))

            self.fields.append(afwTable.Field["F"]("KCDiffimMseResids_%s"%(kType),
                                                   "Mean squared error in diffim : Variance + Bias**2"))

            self.fields.append(afwTable.Field["F"]("KCKernelCentX_%s"%(kType),
                                                   "Centroid in X for this Kernel",
                                                   "pixels"))

            self.fields.append(afwTable.Field["F"]("KCKernelCentY_%s"%(kType),
                                                   "Centroid in Y for this Kernel",
                                                   "pixels"))

            self.fields.append(afwTable.Field["F"]("KCKernelStdX_%s"%(kType),
                                                   "Standard deviation in X for this Kernel",
                                                   "pixels"))

            self.fields.append(afwTable.Field["F"]("KCKernelStdY_%s"%(kType),
                                                   "Standard deviation in Y for this Kernel",
                                                   "pixels"))

            self.fields.append(afwTable.Field["I"]("KernelCandidateId_%s"%(kType),
                                                   "Id for this KernelCandidate"))



        self.fields.append(afwTable.Field["I"]("KCKernelStatus_LOCAL",
                                               "Status of the KernelCandidate"))

        self.fields.append(afwTable.Field["ArrayD"]("KernelCoeffValues_LOCAL",
                                                    "Original basis coefficients",
                                                    nKernelSpatial))

        self.fields.append(afwTable.Field["F"]("BackgroundValue_LOCAL",
                                               "Evaluation of background model at this point"))

        self.fields.append(afwTable.Field["F"]("KCDiffimMseKernel_SPATIAL",
                                               "Mean squared error of spatial kernel estimate"))


    def addToSchema(self, inSourceCatalog):
        """Add the to-be-generated QA keys to the Source schema"""
        schema = inSourceCatalog.getSchema()
        inKeys = []
        fluxKey = inSourceCatalog.getPsfFluxKey()
        centroidKey = inSourceCatalog.getCentroidKey()
        shapeKey = inSourceCatalog.getShapeKey()
        for n in schema.getNames():
            inKeys.append(schema[n].asKey())

        for field in self.fields:
            schema.addField(field)

        outSourceCatalog = afwTable.SourceCatalog(schema)
        for source in inSourceCatalog:
            rec = outSourceCatalog.addNew()
            for k in inKeys:
                if k.getTypeString() == 'Coord':
                    rec.setCoord(source.getCoord())
                else:
                    setter = getattr(rec, "set"+k.getTypeString())
                    getter = getattr(source, "get"+k.getTypeString())
                    setter(k, getter(k))
        outSourceCatalog.definePsfFlux(fluxKey)
        outSourceCatalog.defineCentroid(centroidKey)
        outSourceCatalog.defineShape(shapeKey)
        return outSourceCatalog

    def _calculateStats(self, di, dof=0.):
        """Calculate the core QA statistics on a difference image"""
        mask = di.getMask()
        maskArr = di.getMask().getArray()

        # Create a mask using BAD,SAT,EDGE pixels.  Keep detections
        maskArr &= (mask.getPlaneBitMask("BAD")|mask.getPlaneBitMask("SAT")|mask.getPlaneBitMask("EDGE"))

        # Mask out values based on maskArr
        diArr = ma.array(di.getImage().getArray(), mask=maskArr)
        varArr = ma.array(di.getVariance().getArray(), mask=maskArr)

        # Normalize by sqrt variance, units are in sigma
        diArr /= np.sqrt(varArr)
        mean = diArr.mean()

        # This is the maximum-likelihood extimate of the variance stdev**2
        stdev = diArr.std()
        median = ma.extras.median(diArr)

        # Compute IQR of just un-masked data
        data = ma.getdata(diArr[~diArr.mask])
        iqr = np.percentile(data, 75.) - np.percentile(data, 25.)

        #Calculte chisquare of the residual
        chisq=np.sum(np.power(data, 2.))

        # Mean squared error: variance + bias**2
        # Bias = |data - model| = mean of diffim
        # Variance = |(data - model)**2| = mean of diffim**2
        bias = mean
        variance = np.power(data, 2.).mean()
        mseResids = bias**2 + variance

        # If scipy is not set up, return zero for the stats
        try:
            #In try block because of risk of divide by zero
            rchisq=chisq/(len(data)-1-dof)
            # K-S test on the diffim to a Normal distribution
            import scipy.stats
            D, prob = scipy.stats.kstest(data, 'norm')

            A2, crit, sig = scipy.stats.anderson(data, 'norm')
            # Anderson Darling statistic cand be inf for really non-Gaussian distributions.
            if np.isinf(A2) or np.isnan(A2):
                A2 = 9999.
        except ZeroDivisionError:
            D = 0.
            prob = 0.
            A2 = 0.
            crit = np.zeros(5)
            sig = np.zeros(5)
            rchisq = 0

        return {"mean": mean, "stdev": stdev, "median": median, "iqr": iqr,
                "D": D, "prob": prob, "A2": A2, "crit": crit, "sig": sig,
                "rchisq": rchisq, "mseResids": mseResids}


    def apply(self, candidateList, spatialKernel, spatialBackground, dof=0):
        """Evaluate the QA metrics for all KernelCandidates in the
        candidateList; set the values of the metrics in their
        associated Sources"""
        for kernelCandidate in candidateList:
            source = kernelCandidate.getSource()
            schema = source.schema

            # Calculate ORIG stats (original basis fit)
            if kernelCandidate.getStatus() != afwMath.SpatialCellCandidate.UNKNOWN:
                kType = getattr(diffimLib.KernelCandidateF, "ORIG")
                di = kernelCandidate.getDifferenceImage(kType)
                kernel = kernelCandidate.getKernel(kType)
                kstatus = kernelCandidate.getStatus()
                backgroundValue = kernelCandidate.getBackground(kType)
                kernelValues = kernelCandidate.getKernel(kType).getKernelParameters()
                kernelValues = np.asarray(kernelValues)

                lkim = kernelCandidate.getKernelImage(kType)
                centx, centy = calcCentroid(lkim.getArray())
                stdx, stdy = calcWidth(lkim.getArray(), centx, centy)
                solution = kernelCandidate.getKernelSolution(kType)
                # NOTE
                # What is the difference between kernelValues and solution?

                localResults = self._calculateStats(di, dof=dof)

                metrics = {"KCDiffimMean_LOCAL":localResults["mean"],
                           "KCDiffimMedian_LOCAL":localResults["median"],
                           "KCDiffimIQR_LOCAL":localResults["iqr"],
                           "KCDiffimStDev_LOCAL":localResults["stdev"],
                           "KCDiffimKSD_LOCAL":localResults["D"],
                           "KCDiffimKSProb_LOCAL":localResults["prob"],
                           "KCDiffimADA2_LOCAL":localResults["A2"],
                           "KCDiffimADCrit_LOCAL":localResults["crit"],
                           "KCDiffimADSig_LOCAL":localResults["sig"],
                           "KCDiffimChiSq_LOCAL":localResults["rchisq"],
                           "KCDiffimMseResids_LOCAL":localResults["mseResids"],
                           "KCKernelCentX_LOCAL":centx,
                           "KCKernelCentY_LOCAL":centy,
                           "KCKernelStdX_LOCAL":stdx,
                           "KCKernelStdY_LOCAL":stdy,
                           "KernelCandidateId_LOCAL":kernelCandidate.getId(),
                           "KernelCoeffValues_LOCAL":kernelValues}
                for k in metrics.keys():
                    key = schema[k].asKey()
                    setter = getattr(source, "set"+key.getTypeString())
                    setter(key, metrics[k])
            else:
                try:
                    kType = getattr(diffimLib.KernelCandidateF, "ORIG")
                    lkim = kernelCandidate.getKernelImage(kType)
                except:
                    lkim = None

            # Calculate spatial model evaluated at each position, for
            # all candidates
            skim  = afwImage.ImageD(spatialKernel.getDimensions())
            spatialKernel.computeImage(skim, False, kernelCandidate.getXCenter(),
                                       kernelCandidate.getYCenter())
            centx, centy = calcCentroid(skim.getArray())
            stdx, stdy = calcWidth(skim.getArray(), centx, centy)

            sk = afwMath.FixedKernel(skim)
            sbg = spatialBackground(kernelCandidate.getXCenter(), kernelCandidate.getYCenter())
            di = kernelCandidate.getDifferenceImage(sk, sbg)
            spatialResults = self._calculateStats(di, dof=dof)

            # Kernel mse
            if lkim is not None:
                skim     -= lkim
                bias      = np.mean(skim.getArray())
                variance  = np.mean(np.power(skim.getArray(), 2.))
                mseKernel = bias**2 + variance
            else:
                mseKernel = -99.999

            metrics = {"KCDiffimMean_SPATIAL":spatialResults["mean"],
                       "KCDiffimMedian_SPATIAL":spatialResults["median"],
                       "KCDiffimIQR_SPATIAL":spatialResults["iqr"],
                       "KCDiffimStDev_SPATIAL":spatialResults["stdev"],
                       "KCDiffimKSD_SPATIAL":spatialResults["D"],
                       "KCDiffimKSProb_SPATIAL":spatialResults["prob"],
                       "KCDiffimADA2_SPATIAL":spatialResults["A2"],
                       "KCDiffimADCrit_SPATIAL":spatialResults["crit"],
                       "KCDiffimADSig_SPATIAL":spatialResults["sig"],
                       "KCDiffimChiSq_SPATIAL":spatialResults["rchisq"],
                       "KCDiffimMseResids_SPATIAL":spatialResults["mseResids"],
                       "KCDiffimMseKernel_SPATIAL":mseKernel,
                       "KCKernelCentX_SPATIAL":centx,
                       "KCKernelCentY_SPATIAL":centy,
                       "KCKernelStdX_SPATIAL":stdx,
                       "KCKernelStdY_SPATIAL":stdy,
                       "KernelCandidateId_SPATIAL":kernelCandidate.getId()}
            for k in metrics.keys():
                key = schema[k].asKey()
                setter = getattr(source, "set"+key.getTypeString())
                setter(key, metrics[k])

    def aggregate(self, sourceCatalog, metadata, wcsresids, diaSources=None):
        """Generate aggregate metrics (e.g. total numbers of false
        positives) from all the Sources in the sourceCatalog"""
        for source in sourceCatalog:
            sourceId = source.getId()
            if sourceId in wcsresids.keys():
                #Note that the residuals are not delta RA, delta Dec
                #From the source code "bearing (angle wrt a declination parallel) and distance
                coord, resids = wcsresids[sourceId]
                key = source.schema["RegisterResidualBearing"].asKey()
                setter = getattr(source, "set"+key.getTypeString())
                setter(key, resids[0])
                key = source.schema["RegisterResidualDistance"].asKey()
                setter = getattr(source, "set"+key.getTypeString())
                setter(key, resids[1])
                key = source.schema["RegisterRefPosition"].asKey()
                setter = getattr(source, "set"+key.getTypeString())
                setter(key, afwGeom.Point2D(coord.getRa().asRadians(),
                                            coord.getDec().asRadians()))
        if diaSources:
            metadata.add("NFalsePositivesTotal", len(diaSources))
            nRefMatch = 0
            nSrcMatch = 0
            nunmatched = 0
            for source in diaSources:
                refId = source.get("refMatchId")
                srcId = source.get("srcMatchId")
                if refId > 0:
                    nRefMatch += 1
                if srcId > 0:
                    nSrcMatch += 1
                if refId == 0 and srcId == 0:
                    nunmatched += 1
            metadata.add("NFalsePositivesRefAssociated", nRefMatch)
            metadata.add("NFalsePositivesSrcAssociated", nSrcMatch)
            metadata.add("NFalsePositivesUnassociated", nunmatched)
        for kType in ("LOCAL", "SPATIAL"):
            for sName in ("KCDiffimMean", "KCDiffimMedian", "KCDiffimIQR", "KCDiffimStDev", 
                          "KCDiffimKSProb", "KCDiffimADSig", "KCDiffimChiSq", 
                          "KCDiffimMseResids", "KCDiffimMseKernel"):
                if sName == "KCDiffimMseKernel" and kType == "LOCAL":
                    continue
                kName = "%s_%s" % (sName, kType)
                vals = np.array([s.get(kName) for s in sourceCatalog])
                idx = np.isfinite(vals)
                metadata.add("%s_MEAN" % (kName), np.mean(vals[idx]))
                metadata.add("%s_MEDIAN" % (kName), np.median(vals[idx]))
                metadata.add("%s_STDEV" % (kName), np.std(vals[idx]))

