#!/bin/sh

pwd=`pwd`
# LSST_POLICY_DIR=${pwd}/policy
# export LSST_POLICY_DIR
# echo LSST_POLICY_DIR ${LSST_POLICY_DIR} 

# --------------------------------------------------------- 
# INPUT PARAMETERS
# To run on a single host, keep nodes set equal to 1 
# Increase nodes for a larger parallel execution. 
# For example, for two nodes with 4 cpus we could set nodes=2 
# and nslices=3 (pipeline itself takes one cpu) 
nodes=1
nslices=1
# --------------------------------------------------------- 

# Add 1 to the number of slices to get the universe size 
usize=$(( $nslices + 1 ))

echo "nodes ${nodes}"
echo "nslices ${nslices}"
echo "usize ${usize}"

# MPI commands will be in PATH if mpich2 is in build
echo "Running mpdboot"

mpdboot --totalnum=${nodes} --file=nodelist.scr --verbose

sleep 3s
echo "Running mpdtrace"
mpdtrace -l
sleep 2s

echo "Running mpiexec"

mpiexec -usize ${usize}  -machinefile nodelist.scr -np 1 Pipeline.py

sleep 1s

echo "Running mpdallexit"
mpdallexit

