Mapped surface models (eg: UCLA’s brain cortical surface model) unto a sphere using spherical conformal parametrization with gradient descent. Achieved 45% speedup using conjugate orthogonality constraints.


<img width="447" alt="Screenshot 2024-11-10 at 1 58 25 AM" src="https://github.com/user-attachments/assets/57dd8427-aa62-403d-8aef-4999e55ca8e3">


To build and run the project,

In the main directory, run
 ```
$ rm -r build
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Then execute the following command

```
$ ./main ../brain.obj output
```
