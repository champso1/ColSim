# ColSim

ColSim (Collision Simulation) is a library/simulation tool that allows computation of the cross section of physical processes and generation of events. The program is specifically aimed at processes at the hadron colliders (like the LHC in CERN).


## Building

This library is currently only able to be built on Unix-like systems, i.e. Windows is not supported. If you are on Windows, it is recommended to use WSL.


### Prerequisites

  - CMake >= 3.20
  - gcc (any remotely recent version that supports C++11 should work fine)
  - [LHAPDF](https://www.lhapdf.org/): this library is used to query for PDF information which is required for the calculation of the partonic cross section. The website contains sufficient information to get it downloaded, built, and installed on your system.
  - Gnuplot: used for generating plots of kinematic variables
  - (Optional) Doxygen (and its dependencies): used for generating automatic documentation (not fully finished)

  With the exception of LHAPDF, these can all be installed via your systems package manager without too much issue. As a note: LHAPDF cannot be built on Windows; you must use WSL, or have a Linux/Mac machine (Mac does have some additional installation steps -- see the website).
  
  
### LHAPDF

Once LHAPDF is installed, ensure that all system paths are updated. In particular ensure the 'LD\_LIBRARY\_PATH', 'PATH', and 'PYTHONPATH' environment variables are set correctly, as these are all required for LHAPDF to install PDF sets. Once this is done, simply run

```
lhapdf install <pdf-set>
```

Some defaults in the program are `CT18NNLO` and `cteq6l1`, so these are good to install. To find more, if you want to play around with different sets, go [https://www.lhapdf.org/pdfsets.html](here).

  
### Compiling

As a CMake project, this follows very standard building procedures:

```
mkdir build && cd build
cmake .. -DLHAPDF_ROOT_DIR=<path-to-LHAPDF-install-prefix>
make
make install
```

The default for this project is to install the files in the `install` folder at the top level of the project, not the system paths. You can change it to something else by inserting your install prefix of choice in the third line above, or otherwise leave it blank to install to the default place. If you don't let it install to the default place, you will have to modify the CMakeLists.txt file in the examples to instead point to your chosen installation location.


## Usage

### Input Variables

There are a number of possible input variables that can be specified in a number of ways. Both the hard scattering/cross section calculation and the parton showering parts of the program have different input variables. The variables for the hard/scattering cross section computation are:


- **NumXSIterations**: The number of Monte Carlo iterations to do in the computation of the cross section. This parameter of course carries no physical significance, but higher numbers reduce error and vice versa. The default is one million, but on my machine this already runs in under a second, so higher iterations are very feasible to give very accurate results.
- **ECM**: The center-of-mass energy for the proton-proton collision, or in other words, each proton carries ECM/2. This is set to 14 TeV as this is the current ECM that the LHC in CERN is running collisions at.
- **MinCutoffEnergy**: The energy cutoff imposed during the cros section calculation. There is an absolute cutoff before the entire theory on which these calculations are based breakdown; this is around ~1 GeV. Due to the larger energy scale of the main interaction, we can afford to set this cutoff a bit higher to increase computational efficiency, but this can in principle be set as low or high as one desires, with accuracy impacts.
- **TransformationEnergy**: An intermediate calculational parameter that is usually kept around ~60 GeV, or roughly around the minimum cutoff energy, and used as part of a change-of-variables/transformation to more easily apply the Monte Carlo hit-or-miss method. Do not change this variable below the MinCutoffEnergy or above the square root of the ECM, as this will cause divergences.

The parton showering parameters are as follows:

- **InitialEvolEnergy**: The initial energy that the quark has before undergoing its evolution. This is usually set to the be roughly the energy scale of a generic subprocess of the main proton-proton collision which is roughly ~1 TeV.
- **FixedScale**: This is a Yes/No parameter. The value of the coupling term alpha\_s for the interaction between the quarks and the gluons in principle changes with the energy scale, but not enough in this regime to substantially impact the physics. Of course, letting it vary with the quark's energy scale leads to slightly higher accuracy, but again, the impact is not significant enough to motivate me to keep it one way or the other.
- **EvolutionEnergyCutoff**: As described earlier, this cutoff is due to the fact that at lower energies than this our theory breaks down. Here, since we are directly evolving the quark to these low energies we want to have this cutoff be at the absolute minimum to capture as many emissions as possible. The code will actually error out if the user tries to set this to a value lower than 1 GeV and warns for a value >10 GeV.

These can either be left alone, meaning they take their default values in the program, they can be specified within the program itself via the SETTINGS class; for instance, to change the ECM to 16.7 TeV, one would put:

```cpp
SETTINGS.readString("ECM=16.7");
```

There is also the option for a configuration file. In the `res` directory in the top level of the project is a configuration file titled `config.in`, which contains all of the above variables, their descriptions, and their default values. This can be placed in the same directory as an executable, and by passing its name/path into into the `init()` function of ColSimMain, it will read from their instead. This looks like:

```cpp
colsim.init(ColSimMain::HARD_SCATTERING, "config.in");
```



### Examples

Currently there is only one example located in `examples/basic`, which illustrates the computation of the cross section for the process PP2Zg2ll, and generates a large number events, and generates plots of the kinematic variables. To run this (and any other future examples), navigate to the corresponding directory and invoke the usual CMake commands (minus the installation). You'll be left with an executable whose name matches the example name.


# TODO

Outlined in this section are a handful of some of the next major items on my TODO list to get completed next:

- LHE file interaction.
- Relinquish more plotting functionality to the user.
- In order to do this above point, a cleaner interface to the calculated kinematic variables would be nice.
- Implementation of some more process.
