################################
Alard-Lupton implementation code
################################

This page extends the API documentation of the implementation code.

Solving for the matching kernel starts with the generation of the
basis functions. The number of basis functions and the properties of
the Gaussian bases (e.g. widths in pixels) are determined in
`lsst.ip.diffim.makeKernelBasisList.generateAlardLuptonBasisList`. The
decision logic depends on configuration and psf fwhm values of the
images, some details is given in the function API documentation.

:numref:`figclass` shows the main functional relation among the implementation
classes. This diagram does not detail the abstraction hierarchy of the
classes, they all have more general base classes. Please refer to the
API documentation for details. There are also additional enums and
classes declared to represent solution statuses and to select and
configure numerical algorithms.  ``KernelPca`` and
``KernelPcaVisitor`` are also not detailed here.

.. _figclass:

.. figure:: figures/AL_kernel_solution_classes.svg
    :align: center
    :alt: Kernel solution classes

    Alard-Lupton kernel solution classes and their relationships.

    Arrows point from the class that performs the arrow action on the
    _pointed_ class. All containment or storing relations mean holding
    of *shared_ptrs* of the object instances, there is no exclusive
    ownership.

..
   Figure source on github:  
   lsst-dm/diffimTests/tickets/DM-19443_visualize_coefficients/AL_kernel_solution_classes.drawio

The solution is implemented in `lsst.ip.diffim.PsfMatchTask._solve`
with an iterative loop to select good kernel solutions. In the C++ codebase, the
solution uses the general visitor concept of ``SpatialCellSet`` to
obtain local solutions across all over the image.

In the solution loop, ``BuildSingleKernelVisitor`` visits all kernel
candidates in each ``SpatialCell`` and calls calls their ``build``
method to solve for the coefficents of the basis
kernels. ``KernelCandidate`` owns the knowledge of how to initialize
``KernelSolution``, the solution knows how to solve itself and how to
turn that into an output kernel. The solution is stored in the
``KernelCandidate`` instance itself. Following the solution of each
kernel candidate with the initial basis kernels, a PCA can be
performed to reduce the number of basis functions for the spatial
solution ( ``config.usePcaForSpatialKernel`` ), typically in the case
of the delta-function basis kernels. This needs recalculation of the
local solutions using the new PCA basis kernels and then solving for
the spatial coefficients of the PCA basis. To support PCA,
``KernelCandidate`` can store up to 2 kernel solutions. See C++ API
docs of ``lsst::ip::diffim::KernelCandidate.build``, ``KernelPca`` and
``KernelPcaVisitor`` for more details.

Following the determination of the kernel solutions for each kernel
candidate, the spatial solution is determined by
``BuildSpatialKernelVisitor``. ``BuildSpatialKernelVisitor`` visits
the kernel candidates and collects their static solution and
initializes ``SpatialKernelSolution``. Then ``SpatialKernelSolution``
solves itself.

Both the local _local_ and the whole image _spatial_ solutions are
stored as ``LinearCombinationKernel`` instances. Note, the python API
of ``LinearCombinationKernel`` currently does not support the
evaluation of parameters for an arbitrary x,y position. Workarounds
are to get and evaluate the parameter functions themselves directly or
to compute a kernel image that updates the actual parameter values.
