rm -r decamDirTestNEW
imageDifference2.py calexpDir_b1631 --output decamDirTestNEW \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions >& output_NEW.txt

imageDifference.py calexpDir_b1631 --output decamDirTestOLD \
                   --id visit=289820 ccdnum=11 --templateId visit=288976 \
                   --configfile diffimconfig_OLD.py --clobber-config --clobber-versions >& output_OLD.txt


/bin/cp -fv txt output_NEW.txt decamDirTestNEW/config/deepDiff.py /tmp/
rpl -q makeDiffim. '' /tmp/output_NEW.txt /tmp/deepDiff.py
rpl -q .makeDiffim '' /tmp/output_NEW.txt /tmp/deepDiff.py
rpl -q processDiffim. '' /tmp/output_NEW.txt /tmp/deepDiff.py
rpl -q .processDiffim '' /tmp/output_NEW.txt /tmp/deepDiff.py
diff /tmp/output_NEW.txt output_OLD.txt
diff /tmp/deepDiff.py decamDirTestOLD/config/deepDiff.py

exit

# Metrics dont work (for both versions!)
imageDifference2.py calexpDir_b1631 --output decamDirTestNEW \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
                    --config makeDiffim.doAddMetrics=True

imageDifference.py calexpDir_b1631 --output decamDirTestOLD \
                   --id visit=289820 ccdnum=11 --templateId visit=288976 \
                   --configfile diffimconfig_OLD.py --clobber-config --clobber-versions \
                   --config doAddMetrics=True


# Try Zogy
imageDifference2.py calexpDir_b1631 --output decamDirTestNEWZogy \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile diffimconfig.py --clobber-config --clobber-versions \
                    --config makeDiffim.subtract='zogy' >& output_NEW_ZOGY.txt

imageDifference.py calexpDir_b1631 --output decamDirTestOLDZogy \
                   --id visit=289820 ccdnum=11 --templateId visit=288976 \
                   --configfile diffimconfig_OLD.py --clobber-config --clobber-versions \
                   --config subtract='zogy' >& output_OLD_ZOGY.txt

/bin/cp -fv txt output_NEW_ZOGY.txt decamDirTestNEWZogy/config/deepDiff.py /tmp/
rpl -q makeDiffim. '' /tmp/output_NEW_ZOGY.txt /tmp/deepDiff.py
rpl -q .makeDiffim '' /tmp/output_NEW_ZOGY.txt /tmp/deepDiff.py
rpl -q processDiffim. '' /tmp/output_NEW_ZOGY.txt /tmp/deepDiff.py
rpl -q .processDiffim '' /tmp/output_NEW_ZOGY.txt /tmp/deepDiff.py
diff /tmp/output_NEW_ZOGY.txt output_OLD_ZOGY.txt

exit


# Try individual new command-line tasks. Not sure how they'll work.

makeDiffim.py calexpDir_b1631 --output decamDirTest2 \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile makeDiffimConfig.py --clobber-config --clobber-versions \
                    --config subtract='al'

processDiffim.py decamDirTest2 --output decamDirTest2a \
                    --id visit=289820 ccdnum=11 --templateId visit=288976 \
                    --configfile processDiffimConfig.py --clobber-config --clobber-versions

ls -al decamDirTest2a/deepDiff/v289820/diaSrc-11.fits decamDirTestOLD/deepDiff/v289820//diaSrc-11.fits



# Test with `doUseRegister=True` -- seems to work now that I fixed it
imageDifference2.py calexpDir_b1631 --output decamDirTest_Zogy_noSpatial_doRegister --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile diffimconfig.py --clobber-config --clobber-versions --config makeDiffim.subtract='zogy' --config makeDiffim.doUseRegister=True

