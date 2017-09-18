# Just run the makeDiffim part to get accurate timings of subtraction algorithms/flavors

echo 'AL'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=False >&/tmp/zzz
rm -rf DELETEME
echo 'AL decorr.'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True >&/tmp/zzz
rm -rf DELETEME
echo 'AL decorr. spatial'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True --config doSpatiallyVarying=True >&/tmp/zzz
rm -rf DELETEME

echo 'Zogy'
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy >&/tmp/zzz
rm -rf DELETEME
echo 'Zogy spatial'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy --config doSpatiallyVarying=True >&/tmp/zzz
rm -rf DELETEME
echo 'Zogy im-space'
rpl -q 'inImageSpace=False' 'inImageSpace=True' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy >&/tmp/zzz
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
rm -rf DELETEME
echo 'Zogy im-space spatial'
rpl -q 'inImageSpace=False' 'inImageSpace=True' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=289820 ccdnum=11 --templateId visit=288976 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy --config doSpatiallyVarying=True >&/tmp/zzz
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
rm -rf DELETEME

echo 'AL pre-conv'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=False --config doPreConvolve=True >&/tmp/zzz
rm -rf DELETEME
echo 'AL pre-conv decorr.'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True --config doPreConvolve=True >&/tmp/zzz
rm -rf DELETEME
echo 'AL pre-conv decorr. spatial'
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config doDecorrelation=True --config doPreConvolve=True --config doSpatiallyVarying=True >&/tmp/zzz
rm -rf DELETEME

echo 'Zogy pre-conv'
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy --config doPreConvolve=True >&/tmp/zzz
rm -rf DELETEME
echo 'Zogy pre-conv spatial'
rpl -q 'inImageSpace=True' 'inImageSpace=False' makeDiffimConfig.py
time makeDiffim.py calexpDir_b1631 --output DELETEME --id visit=288976 ccdnum=11 --templateId visit=289820 --configfile makeDiffimConfig.py --clobber-config --clobber-versions --config subtract=zogy --config doPreConvolve=True --config doSpatiallyVarying=True >&/tmp/zzz
rm -rf DELETEME

rm -f /tmp/zzz

