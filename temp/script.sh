cp -fv diffim_config_TEMPLATE.py config.py
echo "config.subtractAlgorithm='AL'" >> config.py
echo "config.doDecorrelation=False" >> config.py
echo "config.detection.thresholdValue=5.5" >> config.py

echo 'A&L'
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=411033 ccdnum=1..3 --templateId visit=410927 \
                                       --output diffim_15A38_g -C config.py --clobber-config \
                                       --clobber-config --no-versions >&diffim_15A38_g.txt

cp -fv diffim_config_TEMPLATE.py config.py
echo "config.subtractAlgorithm='AL'" >> config.py
echo "config.doDecorrelation=True" >> config.py
echo "config.doSpatiallyVarying=False" >> config.py

echo 'A&L decorrelated'
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=411033 ccdnum=1..3 --templateId visit=410927 \
                                       --output diffim_15A38_newDecorr_g -C config.py --clobber-config \
                                       --no-versions >&diffim_15A38_newDecorr_g.txt

cp -fv diffim_config_TEMPLATE.py config.py
echo "config.subtractAlgorithm='AL'" >> config.py
echo "config.doDecorrelation=True" >> config.py
echo "config.doSpatiallyVarying=True" >> config.py

echo 'A&L decorrelated (spatial)'
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=411033 ccdnum=1..3 --templateId visit=410927 \
                                       --output diffim_15A38_newDecorrSpatial_g -C config.py --clobber-config \
                                       --no-versions >&diffim_15A38_newDecorrSpatial_g.txt

cp -fv diffim_config_TEMPLATE.py config.py
echo "config.subtractAlgorithm='ZOGY'" >> config.py
echo "config.doSpatiallyVarying=False" >> config.py

echo 'ZOGY'
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=411033 ccdnum=1..3 --templateId visit=410927 \
                                       --output diffim_15A38_ZOGY_g -C config.py --clobber-config \
                                       --no-versions >&diffim_15A38_ZOGY_g.txt

cp -fv diffim_config_TEMPLATE.py config.py
echo "config.subtractAlgorithm='ZOGY'" >> config.py
echo "config.doSpatiallyVarying=True" >> config.py

echo 'ZOGY (spatial)'
time \
    $PIPE_TASKS_DIR/bin/imageDifference.py processed_15A38 --id visit=411033 ccdnum=1..3 --templateId visit=410927 \
                                       --output diffim_15A38_ZOGYspatial_g -C config.py --clobber-config \
                                       --no-versions >&diffim_15A38_ZOGYspatial_g.txt
