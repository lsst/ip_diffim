.. lsst-task-topic:: lsst.ip.diffim.DiaCatalogSourceSelectorTask

##########################
DiaCatalogSourceSelectorTask
##########################

A task that selects sources for Kernel candidates.

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-description:

Description
==================

A naive star selector based on second moments. Use with caution.

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-debug:


Debug Variables
==================

DiaCatalogSourceSelectorTask has a debug dictionary with the following keys:

display : `bool`
    if True display debug information
displayExposure : `bool`
    if True display exposure
pauseAtEnd `bool`
    if True wait after displaying everything and wait for user input

.. _lsst.ip.diffim.DiaCatalogSourceSelectorTask-examples:

Examples
==================

For example, put something like:

.. code-block:: py

    import lsstDebug
    def DebugInfo(name):
        di = lsstDebug.getInfo(name)  # N.b. lsstDebug.Info(name) would call us recursively
        if name.endswith("diaCatalogSourceSelector"):
            di.display = True

        return di

    lsstDebug.Info = DebugInfo

into your `debug.py` file and run your task with the `--debug` flag.