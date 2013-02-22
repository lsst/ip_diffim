import numpy as np
import matplotlib.pyplot as plt
import sys
from lsst.pipe.tasks.imageDifference import ImageDifferenceTask

# This test uses the kernel sum and spatial model condition number as
# metrics to determine how big the stamps need to be, as a function of
# the kernel size.

# Run like a Task, as in:
# compareKernelSizes.py . --id visit=865833781 raft=2,2 sensor=1,1 --configfile imageDifferenceConfig.py --output=tmplsstdiff 

kSizes = np.arange(13, 33, 2)
gSizes = np.arange(2, 5, 0.25)

kSums  = []
cNums  = []

for kSize in kSizes:
    for gSize in gSizes:
        fpGrowPix = int(gSize * kSize / 2 + 0.5)  # this is a grow radius

        # Specializations for this test
        task_args = sys.argv[1:]
        task_args.append("--config")
        task_args.append("subtract.kernel.active.kernelSize=%d" % (kSize))
        task_args.append("subtract.kernel.active.detectionConfig.fpGrowMin=%d" % (fpGrowPix-1))
        task_args.append("subtract.kernel.active.detectionConfig.fpGrowPix=%d" % (fpGrowPix))
        task_args.append("subtract.kernel.active.detectionConfig.fpGrowMax=%d" % (fpGrowPix+1))
        task_args.append("doDetection=False")
        task_args.append("doMeasurement=False")

        ## Hack around the lack of metadata being returned in cmdLineTask
        tmp = ImageDifferenceTask()
        parsedCmd = tmp._makeArgumentParser().parse_args(config=tmp.ConfigClass(), args=task_args, log=tmp.log, override=tmp.applyOverrides)
        task = ImageDifferenceTask(name = ImageDifferenceTask._DefaultName, config=parsedCmd.config, log=parsedCmd.log)
        results = task.run(parsedCmd.dataRefList[0])

        try:
            kSum = task.subtract.metadata.get("spatialKernelSum")
            cNum = task.subtract.metadata.get("spatialConditionNum")
        except:
            kSums.append(np.inf)
            cNums.append(np.inf)
        else:
            kSums.append(kSum)
            cNums.append(cNum)
        
fig = plt.figure()
data = np.array(kSums).reshape(len(kSizes), len(gSizes)).T
ymin = kSizes.min()
ymax = kSizes.max()
xmin = gSizes.min()
xmax = gSizes.max()
im1 = plt.imshow(data, origin='lower', cmap=plt.cm.jet, extent=[ymin, ymax, xmin, xmax], interpolation="nearest", aspect="auto")
plt.title("Kernel Sum")
plt.xlabel("Kernel Size")
plt.ylabel("Stamp Grow")
plt.colorbar(im1, ax = fig.gca(), orientation="horizontal")

fig = plt.figure()
data = np.array(np.log10(cNums)).reshape(len(kSizes), len(gSizes)).T
ymin = kSizes.min()
ymax = kSizes.max()
xmin = gSizes.min()
xmax = gSizes.max()
im2 = plt.imshow(data, origin='lower', cmap=plt.cm.jet, extent=[ymin, ymax, xmin, xmax], interpolation="nearest", aspect="auto")
plt.title("log10(Condition Number)")
plt.xlabel("Kernel Size")
plt.ylabel("Stamp Grow")
plt.colorbar(im2, ax = fig.gca(), orientation="horizontal")

plt.show()
