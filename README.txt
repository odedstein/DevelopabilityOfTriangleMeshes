This is the code accompanying the paper "Developability of Triangle Meshes" by Oded Stein, Eitan Grinspun, and Keenan Crane.
The paper is to appear in ACM Trans. Graph 2018, vol. 37, number 4.
More info on the paper can be found on its project page, http://www.cs.columbia.edu/cg/developability/


This code is meant to run on macOS 10.13, with Xcode and the Developer tools installed.

This code is built with cmake.
Go to the code directory, and then:

mkdir build
cd build
cmake ..
make

Options can be adjusted with ccmake in the build directory.
Be sure to set CMAKE_BUILD_TYPE to Debug or Release respectively.


The application is run in the build directory by issuing:
./applications/main/main path_to_mesh.obj

There are multiple keypress options to influence the different flows.
A few are listed here:
- TAB: perform one optimization step
- Space: continuously optimize
- up / down arrow keys: increase / decrease step size.
- e: toggles developability energies. 0/1 are the L2/Linf version of the
      covariance energy, modified to reduce spikes. 2/3 are the L2/Linf version
      of the covariance energy sans modification. 4/5 are the L2/Linf version
      of the combinatorial energy.
- m: toggles the minimization strategy. The options are L-BFGS & Overton's
      line search, gradient descent & backtracking linesearch, and gradient
      descent without any linesearch. If a linesearch strategy is active,
      the input from the up / down arrow keys will be overwritten by the
      line search.
- s: in-place 1-4 subdivision of the mesh. Use this, as described in the paper,
      after some time has passed with the flow to get nice, hi-res models.
- n: toggle remeshing of degenerate triangles
- f: cut and flatten mesh. Will write a flattened & cut mesh to
      flattened_out.obj and a mesh with the flattened & cut texture coordinates
      to mesh_with_flattened_texcoord.obj.
      This method is rather sensitive to parameters defined in
      measure_once_cut_twice.cpp at the beginning of the file.
- c: toggles coloring of the mesh
- d: toggles edge coloring by dihedral angle
- w: write the current mesh to out.obj
- v: toggles video recording, this will write a mesh frame*.obj every few steps
- p: toggles the energy plot window
- l: one pass of Laplacian smoothing (tiny timestep)

To run some of our pre-generated examples, run the scripts stored in the directory runs by issuing:
./applications/main/main ../runs/nut.txt


The MATLAB dependency can be removed by altering the respective line in the
CMakeLists.txt from "ON" to "OFF", and changing the preprocessor define on the
top of applications/main/main.cpp from:
#define REQUIRE_MATLAB
to:
#undef REQUIRE_MATLAB
