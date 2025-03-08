# ColSim

ColSim (Collision Simulation) is a library/simulation tool that allows computation of the cross section of physical processes and generation of events. The program is specifically aimed at processes at the hadron colliders (like the LHC in CERN).


## Building

This library is currently only able to be built on Unix-like systems, i.e. Windows is not supported. If you are on Windows, it is recommended to use WSL.


### Prerequisites

  - CMake >= 3.20
  - gcc (any remotely recent version that supports C++11 should work fine)
  - [LHAPDF](https://www.lhapdf.org/): this library is used to query for PDF information which is required for the calculation of the partonic cross section. The website contains sufficient information to get it downloaded, built, and installed on your system. Ensure that your environment variables are set in such a way that the CMake files can be found by this project.
  - (Optional) Doxygen (and its dependencies): used for generating the documentation.
  - (Optional) Gnuplot: used for generating plots of kinematic variables. Not integrated fully yet, so currently this is an optional dependency.
  
  With the exception of LHAPDF, these can all be installed via your systems package manager without too much issue.
  
  
### Compiling

As a CMake project, this follows very standard building procedures:

```
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=<path-to-installation>
cmake --build .
cmake --install .
```

The default for this project is to install the files in the `install` folder at the top level of the project, not the system paths. You can change it to something else by inserting your install prefix of choice in the third line above, or otherwise leave it blank to install to the default place.


## Usage

### Configuration File

The configuration file contains some important information. Namely:

- **ECM**: the center of mass energy in TeV
- **PDFName/PDFMemberNo**: the name and set/member number of the PDF set to use.
- **Process**: The process string. Currently, the default `PP2Zg2ll` (pp initial state, Z-boson/photon intermediate state, dilepton final state) is the only supported one, but in general, the "2"'s separate the initial, intermediate, and final states, and any combination of such states should result in a valid process.
- **NumXSIterations**: The number of Monte Carlo iterations to do in the cross section calculation. I find that 1000000 (the current default) is a good sweet spot for accuracy and runtime.
- **AllowPhotonEmission/AllowGluonEmission**: two yes/no settings for whether to allow the final state particles from the hard scattering process to emit photons or gluons, respectively. Just like with the process string, the current set values are the only allowed options since the only programmed case is gluon emission without photon emission.

> **_NOTE:_**: Ensure that the PDF set is installed on your system. To install the current default (CT18NNLO), use `lhapdf install CT18NNLO` (assuming that you set up LHAPDF correctly). This is a good default, so you can just stick with this one. Others can be found [on their website](https://www.lhapdf.org/pdfsets.html) if you are curious, though.



### Examples

Currently there is only one example located in `examples/basic`, which illustrates the computation of the cross section for the process PP2Zg2ll, and generates 10 random events, printing the most recently generated one to the screen. To run this (and any other future examples), navigate to the corresponding directory and invoke the usual CMake commands (minus the installation). You'll be left with an executable whose name matches the example name. You can play around with configuration file options (only the PDF setname/number, ECM, and number of iterations), or generate more events.


# TODO

Outlined in this section are a handful of some of the next major items on my TODO list to get completed next:

- Finalize LHE file implementation and integration with the main part of the program to allow easy loading of events from/saving events to LHE files. There is a direct consequence to this and it is that it will allow much easier comparison with tools like Pythia8 and MadGraph5.
- Further interaction with the Plotting functionality via Gnuplot to plot histograms of some of the kinematic variables for the process. This also allows more easier comparison with Pythia8.
- Implementation of some more processes.

The issues tab of this repository also contains a few other, more open ended things or goals which are not as pertinent as these, such as specific code organization stuff.
