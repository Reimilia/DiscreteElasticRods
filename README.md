# Discrete Elastic Rods
Final Project of CS395T

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a binary within the folder.

## Run

From within the `build` directory just issue:

    ./ElasticRod_bin

A glfw app should launch displaying a GUI. The code will try to load mesh files from either
the ./meshes or ../meshes folder, so you need to run the binary from either the project root
folder, or a build subdirectory.

## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::viewer::Viewer` (glfw and opengl).

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git

## UI introduction

Mostly it looks similar as what we used in previous assignments.  The binding keys are:

*   'w' for moving forward;
*   's' for moving backward;
*   'a' for moving left;
*   'd' for moving right;
*   ' ' for toggling simulation.

The UI panel is changed, however. I added a test button as usual to justify the program's correctness. The other simulation options will be explained in the report.