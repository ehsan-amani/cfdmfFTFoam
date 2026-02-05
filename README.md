------------------------------------------------------------------------------------
# The Front-Tracking Method (FTM) package for OpenFOAM-9
 Main developer: Dr. Ehsan Amani, Head of the CFDMF group, https://sites.google.com/view/dramani, 2016-2025
 The following CFDMF group students contributed to the old versions of the codes:
 Mahdi Jafari
 Pedram Yousefi
 Mohammad Amin Amini
 Mohammad Bagher Molaei
 Mehran Sharifi

 Developer’s repository link: https://github.com/ehsan-amani/cfdmfFTFoam
 Licensing provisions: GPLv3
 Programming language: C++
 Nature of problem: Two-phase and three-phase flow problems
 Solution method: Tracking interfaces, characterized by surface triangulated meshes, and their interaction
    with a general unstructured Eulerian grid. 
 cfdmfFTFoam has been extended from or inspired by the following codes: 
	1) Ftc3D: Front-Tracking code by Gretar Tryggvason and Sadegh Dabiri, 
	   publicized as a part of the PARIS code (https://doi.org/10.17632/5cb2yrfx7r.1) and
   2) OpenFOAM v9 (www.openfoam.org).

------------------------------------------------------------------------------------
 This package include
 1- "doxygen.tar.xz" : A detailed dixygen manual containing the structure of the library, classes, files, etc.
 2- "src.tar.xz" : The main library
 3- "solvers.tar.xz" : The main solver
 4- "tutorials.tar.xz" : Tutorials, examples, and benchmarks (see the manuscript for the problem definitions)

------------------------------------------------------------------------------------
# Installation: 
 Prerequisite: OpenFOAM Foundation CFD package version 9 (OpenFOAM-9)
 To copy files to correct directories and compile the solver and libraries, open a 
 terminal in the current folder where this readMe file exists, and run:

##use (for one-go installation)
activate of9 environment variables by
of9 #or similar alias in ~/.bashrc to set of9 environment variables
./Allwmake

##or (for step-by-step installation)

###preset
of9 #or similar alias in ~/.bashrc to set of9 environment variables
mkdir -p $WM_PROJECT_USER_DIR/src/lagrangian
mkdir -p $WM_PROJECT_USER_DIR/applications/solvers/multiphase
mkdir -p $WM_PROJECT_USER_DIR/run/multiphase/cfdmf

###extracting
tar -xf src.tar.xz
tar -xf solvers.tar.xz
tar -xf tutorials.tar.xz
tar -xf doxygen.tar.xz

###copying
cp -rp src/frontTracking $WM_PROJECT_USER_DIR/src/lagrangian
cp -rp solvers/* $WM_PROJECT_USER_DIR/applications/solvers/multiphase
cp -rp tutorials/* $WM_PROJECT_USER_DIR/run/multiphase/cfdmf
cp -rp doxygen $WM_PROJECT_USER_DIR/src/lagrangian/frontTracking

###compiling 
cd $WM_PROJECT_USER_DIR/src/lagrangian/frontTracking 
wclean
wmake libso 2>log

cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/cfdmfFTFoam
wclean
wmake 2>log

------------------------------------------------------------------------------------
to run the cases, use ./Allclean and ./Allrun. For instance: 
cd $WM_PROJECT_USER_DIR/run/multiphase/cfdmf/3DDeform/FT3DDeform
./Allclean 
./Allrun 

------------------------------------------------------------------------------------
to run the doxygen manual (do not forget to activate OF9 environment variables by of9):
xdg-open $WM_PROJECT_USER_DIR/src/lagrangian/frontTracking/doxygen/html/index.html

------------------------------------------------------------------------------------
case description: 
3DDeform/FT3DDeform: Front tracking simulation of 3D deformation test
3DDeform/inter3DDeform: VOF simulation of 3D deformation test

RisingBubble/... : Free-rising bubble and rising bubble in a container cases

RisingTaylorBubble/... : Rising Taylor drop in a vertical pipe of circular or cubic cross-sections

ShearDeform/... : Droplet in a Poiseuille flow deformation cases

StagDrop/... : Stagnant droplet cases

stlTocfdmfFTFoamFrontMeshConvertor/ : The binary STL to cfdmfFTFoam front mesh convertor code
