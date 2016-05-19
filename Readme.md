# GraphUtilities library
## Short description

A multi-threaded library for the manipulation of directed graphs. It allows to
perform queries and modifications such as reachability, partitioning, convex
partitioning and coarsening via multiple methods.


## Origin and purpose of the project

This project originated within work done during my Ph.D studies. Its initial
goal was to showcase a novel technique for the indexing of directed acyclic
graphs in order to accelerate reachability queries.

Later on this was pushed further when the need arose to perform convex
partitioning. The convex characteristic implies that the reduced graph of
partitions is also acyclic. 

The next step was to put the convex partitioning to work in the field of the
performance debugging of memory accesses. The studied graphs represent execution
traces of code with memory accesses and reuses and the goal is to compute
alternative schedulings which generate less communications throughout the memory
hierarchy of a processor, i.e less cache misses.


## Using the library
### Install

For the build of the library the CMake tool is used. As such the following steps
should be followed:

1. Create and move to the directory where the build should be performed. In
   order to avoid clobbering the source code I suggest not to use the top-level
   directory of the library code.

2. Use the command-line given below to configure the build. Options are detailed
   further below.
   > `$> cmake <path-to-library-top-level> [OPTIONS]`

3. Build the project with a simple make command
   > `$> make`

4. Optionally create the Doxygen documentation for the library
   > `$> make documentation`

5. Optionally install the library for system-wide use on your computer
   > `$> make install`


Besides the standard CMake options such as `CMAKE_INSTALL_PREFIX`, etc. which
are widely documented online, the GraphUtilities library allows to configure the
following options:

* `-DENABLE_STATISTICS=ON (default: OFF)`
   Builds the library such that it keeps track of query statistics such as
   negative / positive result ratio, search lengths, etc.

* `-DENABLE_BENCHMARKS=ON (default: OFF)`
   Builds the library such that it keeps track of querying and indexation
   execution times. This relies on the PAPI library in order to function.

* `-DENABLE_TLS=OFF (default: ON)`
   Builds the library without using Thread Local Storage. This is used by
   default in order to provide better performances in multi-threaded use-cases.
   I suggest only turning it off if you see really weird behaviour when using
   the library.

* `-DENABLE_COARSEN_STATISTICS (default: OFF)`
   Builds the library such that it keeps track of statistics when performing
   graph coarsening.

* `-DAUTO_THRESHOLD=<int> (default: 10)`
   Amount of queries that should be won by either DFS or BBFS searches in order
   to be selected by automatic search configuration (see the 'Implementation
   details' for more information).

* `-DMAX_THREADS=<int> (default: 2)`
   The maximum number of worker threads that will be used for query processing.
   If the `ENABLE_TLS` options is disabled this may have a significant influence
   on the amount of memory the library uses for storing graphs.


### Integrating the library in your program

The GraphUtilities has an API which is independent from the implementation
through the use of abstract C++ classes.

This API can be observed in the source code in the `src/include` directory, with
the obvious exception of the `src/include/implementation` sub-directory.

For a more clearer oversight I suggest using the Doxygen documentation which can
be accessed by opening `<path-to-build-directory>/doc/public/html/index.html`
in your browser.

A third and last alternative is to build the LaTeX Doxygen documentation in the
`<path-to-build-directory>/doc/public/latex` directory with a `make` command.
The documentation will be build as the `refman.pdf` file.

In order to link to the library you will need to include the
`graph-utilities/graph.hpp` header which can be found in the build directory or
at the install location you specified when configuring CMake.


## Implementation details

Coming soon!
