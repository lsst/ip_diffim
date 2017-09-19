(time imageDifference2.py calexpDir_b1631 --output decamDirTest_AL \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
		    --config makeDiffim.doDecorrelation=False) >& output_AL.txt
cat output_AL.txt

(time imageDifference2.py calexpDir_b1631 --output decamDirTest_ALDec_noSpatial \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
                    --config makeDiffim.doDecorrelation=True) >& output_ALDec_noSpatial.txt
cat output_ALDec_noSpatial.txt

(time imageDifference2.py calexpDir_b1631 --output decamDirTest_ALDec_yesSpatial \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
                    --config makeDiffim.doDecorrelation=True \
		    --config makeDiffim.doSpatiallyVarying=True) >& output_ALDec_yesSpatial.txt
cat output_ALDec_yesSpatial.txt

rpl -q 'inImageSpace=True' 'inImageSpace=False' diffimconfig.py
(time imageDifference2.py calexpDir_b1631 --output decamDirTest_Zogy_noSpatial \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
		    --config makeDiffim.subtract='zogy') >& output_Zogy_noSpatial.txt
cat output_Zogy_noSpatial.txt

rpl -q 'inImageSpace=True' 'inImageSpace=False' diffimconfig.py
(time imageDifference2.py calexpDir_b1631 --output decamDirTest_Zogy_yesSpatial \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
		    --config makeDiffim.subtract='zogy' \
		    --config makeDiffim.doSpatiallyVarying=True) >& output_Zogy_yesSpatial.txt
cat output_Zogy_yesSpatial.txt

rpl -q 'inImageSpace=False' 'inImageSpace=True' diffimconfig.py
(time imageDifference2.py calexpDir_b1631 --output decamDirTest_ZogyImSpace_noSpatial \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
		    --config makeDiffim.subtract='zogy') \
                    >& output_ZogyImSpace_noSpatial.txt
cat output_ZogyImSpace_noSpatial.txt

rpl -q 'inImageSpace=False' 'inImageSpace=True' diffimconfig.py
(time imageDifference2.py calexpDir_b1631 --output decamDirTest_ZogyImSpace_yesSpatial \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
		    --config makeDiffim.subtract='zogy' \
		    --config makeDiffim.doSpatiallyVarying=True) \
                    >& output_ZogyImSpace_yesSpatial.txt
cat output_ZogyImSpace_yesSpatial.txt

rpl -q 'inImageSpace=True' 'inImageSpace=False' diffimconfig.py
