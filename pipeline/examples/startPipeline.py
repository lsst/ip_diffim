from __future__ import with_statement
import sys
import os
import subprocess
import time
import lsst.mwi.utils

def startPipeline(nodeList):
    """Start pipeline execution
    
    Inputs:
    - nodeList: path to node list file

    The pipeline uses one slice just to run the preprocess and postprocess phase;
    all other slices are used to run the processing stages.
    
    Thus to run on a single host specify one host with 2 slices.
    Increase for parallel execution; for example, for two nodes with 4 cpus set nnodes=2 and nslices=3
    
    Warning: requires mpich2 to be setup
    """
    nnodes, nslices = parseNodeList(nodeList)
    usize = nslices + 1 # universe size
    pipelineDir = os.path.dirname(nodeList)
    lsst.mwi.utils.Trace("dps.startPipeline", 3, "nnodes=%s; nslices=%s; usize=%s; pipelineDir=%s" % \
        (nnodes, nslices, usize, pipelineDir))
    
    lsst.mwi.utils.Trace("dps.startPipeline", 3, "Running mpdboot")
    subprocess.call(["mpdboot", "--totalnum=%s" % (nnodes,), "--file=%s" % (nodeList,), "--verbose"])
    
    time.sleep(3)
    lsst.mwi.utils.Trace("dps.startPipeline", 3, "Running mpdtrace")
    subprocess.call(["mpdtrace", "-l"])
    time.sleep(2)
    
    lsst.mwi.utils.Trace("dps.startPipeline", 3, "Running mpiexec")
    subprocess.call(
        ["mpiexec", "-usize", str(usize), "-machinefile", nodeList, "-np", "1", "runPipeline.py"],
        cwd = pipelineDir,
    )
    
    time.sleep(1)    
    lsst.mwi.utils.Trace("dps.startPipeline", 3, "Running mpdallexit")
    subprocess.call(["mpdallexit"])

def parseNodeList(nodeList):
    """Return (nnodes, nslices)"""
    nnodes = 0
    nslices = 0
    with file(nodeList) as nodeFile:
        for line in nodeFile:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue

            try:
                host, slicesStr = line.split(":")
                nslices += int(slicesStr)
            except Exception, e:
                raise RuntimeError("Cannot parse nodeList line %r; error = %s" % (line, e))
            nnodes += 1
    return (nnodes, nslices)
            
            