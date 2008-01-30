import lsst.dps.Stage
import lsst.mwi.data as D
import lsst.mwi.persistence as P
import lsst.fw.Core.fwLib as FW

class VisitMetadataStage(lsst.dps.Stage.Stage):
    def preprocess(self):
        print 'Python VisitMetadataStage preprocess'
        self.activeClipboard = self.inputQueue.getNextDataset()

        inputDP = self.activeClipboard.get("triggerVisitEvent")

        visitDP = D.SupportFactory.createPropertyNode("visit2exposure")
        visitDP.addProperty(D.DataProperty("visitId",
            inputDP.findUnique("visitId").getValueInt()))
        visitDP.addProperty(D.DataProperty.createInt64DataProperty(
            "exposureId", inputDP.findUnique("exposureId").getValueInt() << 1))
        self.activeClipboard.put("visit2exposure", visitDP)

        outputDP = D.SupportFactory.createPropertyNode("rawFPAExposure")

        outputDP.addProperty(D.DataProperty.createInt64DataProperty(
            "rawFPAExposureId",
            inputDP.findUnique("exposureId").getValueInt() << 1))
        outputDP.addProperty(D.DataProperty("ra",
            inputDP.findUnique("FOVRA").getValueDouble()))
        outputDP.addProperty(D.DataProperty("decl",
            inputDP.findUnique("FOVDec").getValueDouble()))
        outputDP.addProperty(D.DataProperty.createFloatDataProperty("equinox",
            inputDP.findUnique("equinox").getValueDouble()))
        outputDP.addProperty(D.DataProperty("filterId",
            self._lookup(inputDP.findUnique("filterName").getValueString())))
        outputDP.addProperty(D.DataProperty.createDateTimeDataProperty(
            "dateObs",
            P.DateTime(inputDP.findUnique("visitTime").getValueDouble())
            ))
        outputDP.addProperty(D.DataProperty.createFloatDataProperty(
            "expTime", inputDP.findUnique("exposureTime").getValueDouble()))
        outputDP.addProperty(D.DataProperty.createFloatDataProperty(
            "airmass", inputDP.findUnique("airMass").getValueDouble()))

        self.activeClipboard.put("rawFPAExposure", outputDP)

	# Rely on postprocess() to put the activeClipboard on the outputQueue.
        
    def _lookup(self, filterName):
        dbLocation = P.LogicalLocation( \
                "mysql://lsst10.ncsa.uiuc.edu:3306/test")
        filterDB = FW.Filter(dbLocation, filterName)
        filterId = filterDB.getId()
        return filterId
