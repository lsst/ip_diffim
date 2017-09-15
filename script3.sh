# Just run the makeDiffim part to get accurate timings of subtraction algorithms/flavors

time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=False >&/tmp/zzz
rm -rf DELETEME
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True >&/tmp/zzz
rm -rf DELETEME
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True --config doSpatiallyVarying=True >&/tmp/zzz
rm -rf DELETEME
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy >&/tmp/zzz
rm -rf DELETEME
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy --config doSpatiallyVarying=True >&/tmp/zzz
rm -rf DELETEME
rpl -q 'inImageSpace=False' 'inImageSpace=True' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy >&/tmp/zzz
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
rm -rf DELETEME

time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=False --config doPreConvolve=True >&/tmp/zzz
rm -rf DELETEME
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True --config doPreConvolve=True >&/tmp/zzz
rm -rf DELETEME

rm -f /tmp/zzz

