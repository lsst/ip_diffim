################################
Alard-Lupton implementation code
################################

This page focuses on key implementation details of the Alard-Lupton
(AL) image differencing in the LSST pipeline.

Introduction
------------

In image differencing, our goal is to detect variable brightness
features between two images by removing their static brightness parts
by subtraction. Before performing pixel by pixel subtraction, the
images should be transformed to match their PSFs. A suitable
transformation is convolution of one of the images by a small image, a
convolution kernel. The optimal convolution kernel is called the
*matching kernel*. Given the images :math:`R` and :math:`I`, we want
to minimise the following expression in the least squares sense by
finding :math:`K`:

.. math::
   R\otimes K - I

In the AL approach the matching kernel is searched in the form of a
linear combination of given basis functions:

.. math::
   K = \sum_i a_i K_i

As the PSFs are not uniform all over the images, the coefficients
:math:`a_i` are considered to be functions of image coordinates. In
the AL algorithm, first local kernel solutions (spatially not varying
values of :math:`a_i`-s) are determined around selected sources. Then
polynomials are fitted to the local sets of :math:`a_i` solutions
to obtain the spatially dependent coefficients for the image pair.

Basis functions
---------------

Solving for the matching kernel starts with the generation of the
basis functions. The number of basis functions and the properties of
the Gaussian bases (e.g. widths in pixels) are determined in
`lsst.ip.diffim.makeKernelBasisList.generateAlardLuptonBasisList`. The
decision logic depends on configuration and psf fwhm values of the
images, some details is given in the function API documentation.

.. TODO More details about the basis function generation

Implementation of the AL solution
---------------------------------

In the LSST pipeline, the following main steps can be distinguished:

* A subset of detected or external catalog sources are selected on a
  spatial grid in the image. *Kernel candidates* (``KernelCandidate``
  in :numref:`figclass`) are created to perform the determination of
  the local matching kernel around each source. The entry point for
  the kernel optimization is at `lsst.ip.diffim.PsfMatchTask._solve()`
  in the Python codebase.

* For each kernel candidate, local (spatially non-varying)
  coefficients are determined of the basis kernels that minimizes the
  image difference for the image stamp around this source.

* If principal component analysis (PCA) is selected, new (fewer) basis
  functions are calculated from the local solutions, then the kernel
  candidates are solved again with the new set of basis functions.

* For the whole image, spatially varying coefficients are fitted onto
  the local coefficient solutions. This is an iterative process
  rejecting outlying local solutions. The spatial variation is either
  generic 2nd order polynomials, or Chebyshev polynomials of the first
  kind (``config.spatialModelType``).

:numref:`figclass` shows the main functional relation among the
implementation classes. These classes are implemented in the C++
codebase. Each arrow represents a functional relationship, pointing
from the subject to its target (object) class. The targets are
created, contained or accessed (used) by the subject classes. Numbers
show the possible multiplicities of the target instances in each
action. This diagram show the functionally most important components
only. It does not detail the abstraction hierarchy of the classes, all
classes shown have more general base classes. There are also
additional enums and classes declared to represent internal statuses
and to select and configure numerical algorithms. Classes for the
statistical analysis of the local and spatial solutions and for
supporting PCA (``KernelSumVisitor``, ``AssessSpatialKernelVisitor``,
``KernelPca``, ``KernelPcaVisitor``) are also not detailed
here. Please refer to the `C++ API
<http://doxygen.lsst.codes/stack/doxygen/x_masterDoxyDoc/index.html>`_
documentation for these details.

.. _figclass:

.. figure:: figures/AL_kernel_solution_classes.svg
    :align: center
    :alt: Kernel solution classes

    Main Alard-Lupton kernel solution classes in the C++ codebase and
    their relationships.

    Arrows point from the class that performs the arrow action on the
    *pointed* class. Numbers mark target multiplicity, "0.." marks
    zero or more, "1.." marks one or more target instances. All
    containment or storing relations mean holding of *shared_ptrs* of
    the object instances, there is no exclusive ownership.  Figure
    source :download:`here
    <figures/AL_kernel_solution_classes.drawio>`.

In the C++ codebase, the solution steps are implemented using the
*visitor concept*. The visitor classes visit each ``KernelCandidate``
instance, perform their operation on the kernel candidate itself by
changing the state of the *visited instance* while also updating their
own *visitor* state with data about the whole visiting process like
the number of good or bad numerical solutions.

``BuildSingleKernelVisitor`` visits all kernel candidates in each
``SpatialCell`` and calls their ``build()`` method to solve for
the coefficients of the basis kernels. ``KernelCandidate`` owns the
knowledge of how to initialize ``KernelSolution``, the solution knows
how to solve itself and how to turn that into an output kernel. The
solution is stored in the ``KernelCandidate`` instance
itself. Following the solution of each kernel candidate with the
initial basis kernels, a principal component analysis (PCA) can be
performed to reduce the number of basis functions for the spatial
solution ( ``config.usePcaForSpatialKernel`` ), typically in the case
of the delta-function basis kernels. This needs the calculation of the
PCA basis of the initial local solutions, recalculation of the local
solutions using the new PCA basis kernels and then solving for the
spatial coefficients of the PCA basis. To support PCA,
``KernelCandidate`` can store one *original* and one *PCA* kernel
solution. See the C++ API docs of ``KernelCandidate.build()``,
``KernelPca`` and ``KernelPcaVisitor`` for more details.

Following the determination of the kernel solutions for each kernel
candidate, the spatial solution is determined by
``BuildSpatialKernelVisitor``. ``BuildSpatialKernelVisitor`` visits
the kernel candidates and collects their local solution and
initializes ``SpatialKernelSolution``. Then ``SpatialKernelSolution``
solves itself.

Both the *local* and the whole image *spatial* solutions are stored as
``LinearCombinationKernel`` instances. Note, the python API of
``LinearCombinationKernel`` currently does not support the evaluation
of parameters for an arbitrary x,y position. Workarounds are to get
and evaluate the parameter functions themselves directly or to compute
a kernel image that updates the last parameter values retrievable from
the instance.
