![alt text](VMTutorial.png)

Lorentz Centre workshop Computational Advances in Active Matter (11 - 15 December 2023)
The Erwin Schrödinger International Institute for Mathematics and Physics summer school on Non-equilibrium Processes in Physics and Biology
(19 - 30 August 2024)

# Tutorial on Vertex Model

A showcase implementation of the vertex model for tissue mechanics.

## Project description

This tutorial provides a simple C++ implementation of the basic vertex model for tissue mechanics. 
Python interface is provided using the [pybind11](https://github.com/pybind/pybind11) library.

The aim of this tutorial is to be clear and pedagogical even it if it comes at some performance cost. 
This tutorial is not intended for research, but could be used to develop research-quality tools.

## Requirements

The package is tested on Linux and Mac OSX.

You will need:

- boost 
- pyvista 

Depending on your local Python installation, both can be installed through conda or pip.

## Installation

Clone the code repository into the VMTutorial directory.

From the VMTutorial directory type:

```
pip install .
```

This should build and install the package into your local site-packages directory. This may take a while since VTK library is cloned and
compiled locally. This is due to common problems with system installations of the VTK library. If, however, a working installation of VTK 
is available, one can add

```
export VTK_DIR=/path/to/vtk/cmake
```

before calling 'pip'

For example, if VTK 9.2 is installed in $HOME/software/VTK/9.2, one would set VTK_DIR=$HOME/software/VTK/9.2/lib/cmake/vtk-9.2

This should significantly speed up the build process.

## Structure

- examples - contains several examples on how to run a simulation
- VMToolkit - the source code for the vertex model simulation and data analysis
- config_builder - several tools for building initial configurations

## Running 

If the installation was successful, one can type (e.g., is an IPython console):

```
from VMToolkit.VM import *
from VMToolkit.VMAnalysis import *
```
