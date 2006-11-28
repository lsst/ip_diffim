#!/bin/sh

# setup all the nodes for mpi testing of RAM/local disk use

rm -rf /dev/shm/MixTMPFS/* /dev/shm/TMPFSExImg/* /dev/shm/UseTMPFS/*
ssh ds29 rm -rf /dev/shm/MixTMPFS/* /dev/shm/TMPFSExImg/* /dev/shm/UseTMPFS/*
ssh ds28 rm -rf /dev/shm/MixTMPFS/* /dev/shm/TMPFSExImg/* /dev/shm/UseTMPFS/*
ssh ds27 rm -rf /dev/shm/MixTMPFS/* /dev/shm/TMPFSExImg/* /dev/shm/UseTMPFS/*
ssh ds26 rm -rf /dev/shm/MixTMPFS/* /dev/shm/TMPFSExImg/* /dev/shm/UseTMPFS/*

/lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpiRamSetup.sh
ssh ds29 /lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpiRamSetup.sh
ssh ds28 /lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpiRamSetup.sh
ssh ds27 /lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpiRamSetup.sh
ssh ds26 /lsst_ibrix/rallsman/LSSTProto/Pipeline/WCS/wcsTest/mpi/mpiRamSetup.sh

