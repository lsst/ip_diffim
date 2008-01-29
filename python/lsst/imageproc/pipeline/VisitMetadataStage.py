import lsst.dps.Stage
import lsst.mwi.data as D
import lsst.mwi.persistence as P

class VisitMetadataStage(lsst.dps.Stage.Stage):
    def preprocess(self):
        print 'Python VisitMetadataStage preprocess'
        activeClipboard = self.inputQueue.getNextDataset()

        inputDPName = _policy.getString("InputDataProperty")
        outputDPName = _policy.getString("OutputDataProperty")

        inputDP = activeClipboard.get("triggerVisitEvent")
        outputDP = D.SupportFactory.createPropertyNode("visitMetadata")

        outputDP.addProperty(D.DataProperty("rawFPAExposureId",
            inputDP.findUnique("exposureId").getValueInt()))
        outputDP.addProperty(D.DataProperty("ra",
            inputDP.findUnique("FOVRA").getValueDouble()))
        outputDP.addProperty(D.DataProperty("decl",
            inputDP.findUnique("FOVDec").getValueDouble()))
        outputDP.addProperty(D.createFloatDataProperty("equinox",
            inputDP.findUnique("equinox").getValueDouble()))
        # \todo Add lookup function
        outputDP.addProperty(D.DataProperty("filterId",
            lookup(inputDP.findUnique("filterName").getValueString())))
        outputDP.addProperty(D.DataProperty.createDateTimeDataProperty( \
            "dateObs", \
            P.DateTime(inputDP.findUnique("visitTime").getValueDouble())
            ))
        outputDP.addProperty(D.DataProperty.createFloatDataProperty(
            "expTime", inputDP.findUnique("exposureTime").getValueDouble()))
        outputDP.addProperty(D.DataProperty.createFloatDataProperty(
            "airmass", inputDP.findUnique("airmass").getValueDouble()))

        activeClipboard.put("visitMetadata", outputDP)
        
        self.outputQueue.addDataset(activeClipboard)

# OutputStage Policy
#
# Formatter: {
#   DataProperty: {
#     visitMetadata: {
#       TableName: "Raw_FPA_Exposure"
#       KeyList: "rawFPAExposureId" "ra" "decl" "equinox"
#       KeyList: "filterId" "dateObs" "expTime" "airmass"
#     }
#   }
# }
