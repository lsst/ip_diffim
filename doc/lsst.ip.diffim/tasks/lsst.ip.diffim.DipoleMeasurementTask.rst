.. lsst-task-topic:: lsst.ip.diffim.DipoleMeasurementTask

##########################
DipoleMeasurementTask
##########################
The list of badFlags will be used to make a list of keys to check for measurement flags on.  By
default the centroid keys are added to this list

.. _lsst.ip.diffim.DipoleMeasurementTask-description:

Description
==================

This class provides a default configuration for running Source measurement on image differences.

.. code-block:: py

    class DipoleMeasurementConfig(SingleFrameMeasurementConfig):
        "Measurement of detected diaSources as dipoles"
        def setDefaults(self):
            SingleFrameMeasurementConfig.setDefaults(self)
            self.plugins = ["base_CircularApertureFlux",
                            "base_PixelFlags",
                            "base_SkyCoord",
                            "base_PsfFlux",
                            "ip_diffim_NaiveDipoleCentroid",
                            "ip_diffim_NaiveDipoleFlux",
                            "ip_diffim_PsfDipoleFlux",
                            "ip_diffim_ClassificationDipole",
                            ]
            self.slots.calibFlux = None
            self.slots.modelFlux = None
            self.slots.instFlux = None
            self.slots.shape = None
            self.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
            self.doReplaceWithNoise = False

These plugins enabled by default allow the user to test the hypothesis that the Source is a dipole.
This includes a set of measurements derived from intermediate base classes
DipoleCentroidAlgorithm and DipoleFluxAlgorithm.
Their respective algorithm control classes are defined in
DipoleCentroidControl and DipoleFluxControl.
Each centroid and flux measurement will have _neg (negative)
and _pos (positive lobe) fields.

The first set of measurements uses a "naive" algorithm
for centroid and flux measurements, implemented in
NaiveDipoleCentroidControl and NaiveDipoleFluxControl.
The algorithm uses a naive 3x3 weighted moment around
the nominal centroids of each peak in the Source Footprint.  These algorithms fill the table fields
ip_diffim_NaiveDipoleCentroid* and ip_diffim_NaiveDipoleFlux*

The second set of measurements undertakes a joint-Psf model on the negative
and positive lobe simultaneously. This fit simultaneously solves for the negative and positive
lobe centroids and fluxes using non-linear least squares minimization.
The fields are stored in table elements ip_diffim_PsfDipoleFlux*.

Because this Task is just a config for SingleFrameMeasurementTask, the same result may be achieved by
manually editing the config and running SingleFrameMeasurementTask. For example:

.. code-block:: py

    config = SingleFrameMeasurementConfig()
    config.plugins.names = ["base_PsfFlux",
                            "ip_diffim_PsfDipoleFlux",
                            "ip_diffim_NaiveDipoleFlux",
                            "ip_diffim_NaiveDipoleCentroid",
                            "ip_diffim_ClassificationDipole",
                            "base_CircularApertureFlux",
                            "base_SkyCoord"]

    config.slots.calibFlux = None
    config.slots.modelFlux = None
    config.slots.instFlux = None
    config.slots.shape = None
    config.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
    config.doReplaceWithNoise = False

    schema = afwTable.SourceTable.makeMinimalSchema()
    task = SingleFrameMeasurementTask(schema, config=config)-



The ``pipetask`` command line interface supports a
flag ``--debug`` to import ref:  supports a flag `--debug`` to import debug.py from your PYTHONPATH; see :ref:`lsstDebug`
for more about debug.py files.

The available variables in DipoleMeasurementTask include:

.. code-block:: py

    import sys
    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)
        if name == "lsst.ip.diffim.dipoleMeasurement":
            di.display = True                 # enable debug output
            di.maskTransparency = 90          # display mask transparency
            di.displayDiaSources = True       # show exposure with dipole results
        return di
    lsstDebug.Info = DebugInfo
    lsstDebug.frame = 1

    config.slots.calibFlux = None
    config.slots.modelFlux = None
    config.slots.gaussianFlux = None
    config.slots.shape = None
    config.slots.centroid = "ip_diffim_NaiveDipoleCentroid"
    config.doReplaceWithNoise = False

Start the processing by parsing the command line, where the user has the option of
enabling debugging output and/or sending their own image for demonstration
(in case they have not downloaded the afwdata package).

.. code-block:: py

    if __name__ == "__main__":
        import argparse
        parser = argparse.ArgumentParser(
            description="Demonstrate the use of SourceDetectionTask and DipoleMeasurementTask")
        parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
        parser.add_argument("--image", "-i", help="User defined image", default=None)
        args = parser.parse_args()
        if args.debug:
            try:
                import debug
                debug.lsstDebug.frame = 2
            except ImportError as e:
                print(e, file=sys.stderr)
        run(args)

The processing occurs in the run function.  We first extract an exposure from disk or afwdata, displaying
it if requested:

.. code-block:: py

    def run(args):
        exposure = loadData(args.image)
        if args.debug:
            afwDisplay.Display(frame=1).mtv(exposure)

Create a default source schema that we will append fields to as we add more algorithms:

.. code-block:: py

    schema = afwTable.SourceTable.makeMinimalSchema()

Create the detection and measurement Tasks, with some minor tweaking of their configs:

.. code-block:: py

        # Create the detection task
    config = SourceDetectionTask.ConfigClass()
    config.thresholdPolarity = "both"
    config.background.isNanSafe = True
    config.thresholdValue = 3
    detectionTask = SourceDetectionTask(config=config, schema=schema)
    # And the measurement Task
    config = DipoleMeasurementTask.ConfigClass()
    config.plugins.names.remove('base_SkyCoord')
    algMetadata = dafBase.PropertyList()
    measureTask = DipoleMeasurementTask(schema, algMetadata, config=config)

Having fully initialised the schema, we create a Source table from it:

.. code-block:: py

    # Create the output table
    tab = afwTable.SourceTable.make(schema)

Run detection:

.. code-block:: py

    # Process the data
    results = detectionTask.run(tab, exposure)

Because we are looking for dipoles, we need to merge the positive and negative detections:

.. code-block:: py

    # Merge the positive and negative sources
    fpSet = results.fpSets.positive
    growFootprint = 2
    fpSet.merge(results.fpSets.negative, growFootprint, growFootprint, False)
    diaSources = afwTable.SourceCatalog(tab)
    fpSet.makeSources(diaSources)
    print("Merged %s Sources into %d diaSources (from %d +ve, %d -ve)" % (len(results.sources),
                                                                      len(diaSources),
                                                                      results.fpSets.numPos,
                                                                      results.fpSets.numNeg))

Finally, perform measurement (both standard and dipole-specialized) on the merged sources:

.. code-block:: py

    measureTask.run(diaSources, exposure)

Optionally display debugging information:

.. code-block:: py

    # Display dipoles if debug enabled
    if args.debug:
        dpa = DipoleAnalysis()
        dpa.displayDipoles(exposure, diaSources)
