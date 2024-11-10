Mapped surface models (eg: UCLA’s brain cortical surface model) unto a sphere using spherical conformal parametrization with gradient descent. Achieved 45% speedup using conjugate orthogonality constraints.

To build and run the project,

In the main directory, run

$ rm -r build
$ mkdir build
$ cd build
$ cmake ..
$ make

Then execute the following command

$ ./main ../brain.obj output
