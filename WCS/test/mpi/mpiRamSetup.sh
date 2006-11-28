#!/bin/sh
mkdir -p /dev/shm/tests/WCS/data
mkdir -p /dev/shm/etc
cp /lsst_ibrix/LSSTProto/tests/WCS/data/642538p.fits /dev/shm/tests/WCS/data/
cp /lsst_ibrix/LSSTProto/etc/*.conf /dev/shm/etc/
mkdir -p /dev/shm/MixTMPFS
mkdir -p /dev/shm/TMPFSExImg
mkdir -p /dev/shm/UseTMPFS

cp /lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpi*.py /dev/shm/UseTMPFS/
cp /lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpi*.py /dev/shm/TMPFSExImg
