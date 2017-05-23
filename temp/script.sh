set CCDNUMS = '1..3'
set VISIT = '411033'
set TEMPLATE = '410927'

cp -fv diffim_config_TEMPLATE.py config1.py
echo "config.subtractAlgorithm='AL'" >> config1.py
echo "config.doDecorrelation=False" >> config1.py
echo "config.detection.thresholdValue=5.5" >> config1.py

echo 'A&L'
grep -v '#' config1.py
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=${VISIT} ccdnum=${CCDNUMS} --templateId visit=${TEMPLATE} \
                                       --output diffim_15A38_g -C config1.py --clobber-config \
                                       --clobber-config --no-versions >&diffim_15A38_g.txt  &

cp -fv diffim_config_TEMPLATE.py config2.py
echo "config.subtractAlgorithm='AL'" >> config2.py
echo "config.doDecorrelation=True" >> config2.py
echo "config.doSpatiallyVarying=False" >> config2.py

echo 'A&L decorrelated'
grep -v '#' config2.py
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=${VISIT} ccdnum=${CCDNUMS} --templateId visit=${TEMPLATE} \
                                       --output diffim_15A38_newDecorr_g -C config2.py --clobber-config \
                                       --no-versions >&diffim_15A38_newDecorr_g.txt  &

cp -fv diffim_config_TEMPLATE.py config3.py
echo "config.subtractAlgorithm='AL'" >> config3.py
echo "config.doDecorrelation=True" >> config3.py
echo "config.doSpatiallyVarying=True" >> config3.py

echo 'A&L decorrelated (spatial)'
grep -v '#' config3.py
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=${VISIT} ccdnum=${CCDNUMS} --templateId visit=${TEMPLATE} \
                                       --output diffim_15A38_newDecorrSpatial_g -C config3.py --clobber-config \
                                       --no-versions >&diffim_15A38_newDecorrSpatial_g.txt

cp -fv diffim_config_TEMPLATE.py config4.py
echo "config.subtractAlgorithm='ZOGY'" >> config4.py
echo "config.doSpatiallyVarying=False" >> config4.py

echo 'ZOGY'
grep -v '#' config4.py
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=${VISIT} ccdnum=${CCDNUMS} --templateId visit=${TEMPLATE} \
                                       --output diffim_15A38_ZOGY_g -C config4.py --clobber-config \
                                       --no-versions >&diffim_15A38_ZOGY_g.txt  &

cp -fv diffim_config_TEMPLATE.py config5.py
echo "config.subtractAlgorithm='ZOGY'" >> config5.py
echo "config.doSpatiallyVarying=True" >> config5.py

echo 'ZOGY (spatial)'
grep -v '#' config5.py
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=${VISIT} ccdnum=${CCDNUMS} --templateId visit=${TEMPLATE} \
                                       --output diffim_15A38_ZOGYspatial_g -C config5.py --clobber-config \
                                       --no-versions >&diffim_15A38_ZOGYspatial_g.txt

rm -v config[12345].py
