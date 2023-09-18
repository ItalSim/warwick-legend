# warwick-legend
Legend muon background simulations

First attempt at creating a LEGEND Monte-Carlo simulation application,
suitable for large production runs estimating Germanium-77 production in
Germanium-76 crystals in different underground laboratories.

Should allow for two alternative experimental setups.

## Requirements

Multi-threading operation and MPI capability

ROOT ntuple output

Biasing of neutron production by muons (physics bias),
neutron population inside cryostant (geometry bias) and
neutron capture reaction (Ge-77 production from n,gamma reaction)

Operate with command line input and Geant4 macro.

Macro commands for change of geometry, 
primary vertex generator and cross section bias factor (one each, neutron and muon)

CLI for number of threads, macro file, output file name (for production runs)

User Limits for run time optimization

## How to Build/Develop
The project has the following requirements:

- Linux/macOS system (only CentOS7, Catalina tested at present)
- C++17 compatible compiler (GCC9, Xcode 11 tested at present)
- CMake 3.12 or newer
- Geant4 10.6 built with multithreading and gdml support
  - Qt or OpenGL/X11 driver support if you want to visualize the geometry/tracks/hits

It may be compiled using:

```console
$ mkdir build && cd build
$ cmake ..
$ make
```

Testing is enabled by default via the `BUILD_TESTING` CMake argument and may be run after
building by running

```
$ ctest
```

Run `ctest --help` for options to select tests and/or run with increased verbosity, e.g.

```
$ ctest -VV
```

to get verbose (including stdout) output from the setup and execution of the tests.

The resulting `warwick-legend` application may be run without arguments to start an interactive
session. Otherwise run `warwick-legend --help` to see a list of options for batch mode running.


Support files for `[clang-format](https://clang.llvm.org/docs/ClangFormat.html)` and `[clang-tidy](https://clang.llvm.org/extra/clang-tidy/)` are provided to help automate formatting (following the
Geant4 style) and static analysis respectively. No explicit setup is needed
for autoformatting other than the relevant integration settings of your
favoured editor/IDE ([vim](http://clang.llvm.org/docs/ClangFormat.html#vim-integration), [emacs](http://clang.llvm.org/docs/ClangFormat.html#emacs-integration), [vscode](https://code.visualstudio.com/docs/cpp/cpp-ide#_code-formatting)). To enable static
analysis, ensure you have `clang-tidy` installed and run the build as:

```console
$ mkdir build-tidy && cd build-tidy
$ cmake -DCMAKE_CXX_CLANG_TIDY="/path/to/clang-tidy" ..
$ make
```

The `.clang-tidy` file supplied in this project will be used, and suggestions
for fixes will be emitted whilst building. At present, the `-fix` option is
automatically apply the suggested change is not used to leave the decision
up to the developer. The set of fixes applied are:

- `readability-*`
- `modernize-*`
- `performance-*`

For a full listing of the wildcards, see [the `clang-tidy` documentation](https://clang.llvm.org/extra/clang-tidy/checks/list.html).

## Code Details
### Cross section bias reference
M.H. Mendenhall and R.A. Weller, NIM A667 (2012) 38-43

### Ntuple output columns
- Hit data, one row per event
  - Edep
  - Time
  - Weight
  - Hit x location
  - Hit y location
  - Hit z location
- Trajectory data, one row per event
  - PDG code
  - N entries in position containers
  - Vertex logical volume name code, see name map
  - Vertex x location
  - Vertex y location
  - Vertex z location
- Trajectory data track step points, N rows, see N entries column above
  - x position
  - y position
  - z position

### Vertex Name Map
Volume definitions in detector construction.
- lookup["Cavern_log"]   = 0;
- lookup["Hall_log"]     = 1;
- lookup["Tank_log"]     = 2;
- lookup["Water_log"]    = 3;
- lookup["Cout_log"]     = 4;
- lookup["Cvac_log"]     = 5;
- lookup["Cinn_log"]     = 6;
- lookup["Lar_log"]      = 7;
- lookup["Lid_log"]      = 8;
- lookup["Bot_log"]      = 9;
- lookup["Copper_log"]   = 10;
- lookup["ULar_log"]     = 11;
- lookup["Ge_log"]       = 12;
- lookup["Pu_log"]       = 13;
- lookup["Membrane_log"] = 14;



# Changes after the fork from the original

A lot has been added in my version of the fork. 
An output for the Ge77m veto has been added, more options to adjust the geometry via macros and several new primary generators.

## Overview over the Macros

### Runaction Macro
Macros regarding the output of the simulation 
```
/WLGD/runaction/
  - WriteOutNeutronProductionInfo
  - WriteOutGeneralNeutronInfo
  - getIndividualGeDepositionInfo
  - getIndividualGdDepositionInfo
```
### Event Macro
Macro to adjust the condition to save all events (1) or just the ones with Ge77 production (0) 
```
/WLGD/event/
  - saveAllEvents (no: [0], yes: 1)
```
### Generator Macro
Macros to controll the primary generator 
```
/WLGD/generator/
  - depth
  - setMUSUNFile (path to file)
  - setMUSUNDirectory (full path to directory containing MUSUN files)
  - setGenerator (options: "MeiAndHume", "Musun", "Ge77m", "Ge77andGe77m", "ModeratorNeutrons", "ExternalNeutrons")
```

More information on SetMUSUNDirectory can be found in the OpenMUSUNDirectory method of src/WLGDPrimaryGeneratorAction.cc

### Detector Macro
Macros to controll the detector geometry
```
/WLGD/detector/
  - setPositionOfDetectors
  - setGeometry
  - exportGeometry
  - XeConc
  - He3Conc
  - Cryostat_Radius_Outer
  - Cryostat_Height
  - Without_Cupper_Tubes
  - With_Gd_Water
  - With_NeutronModerators (options: 0: [no moderators], 1: around re-entrance tubes, 2: in turbine mode, 3: in large hollow tube mode, 5: n-sided hollow polyhedron shield)
  - Which_Material (options: [BoratedPE], PolyEthylene, PMMA plus additional options for doped PMMA can be found in src/WLGDDetectorConstruction.cc)
- TurbineAndTube
    - Radius
    - Width
    - Height
    - zPosition
    - NPanels
    - NSides (for polyhedral tube only)
```
### Bias Macro
Macros to adjust the bias of the cross-sections
```
/WLGD/bias/
  - setNeutronBias
  - setMuonBias
  - setNeutronYieldBias
```
### Step Macro
Macros to adjust whether additional output (additional to the Ge77 production) is recorded in the first place
```
/WLGD/step/
  - getDepositionInfo (multiplicity and energy deposition in the detectors)
  - getIndividualDepositionInfo (energy depositions in the whole cryostat)
  - AllowForLongTimeEmissionReadout (allow for energy depositions >1s after muon crossing to be recorded)
```
## Example for Ge77 production by Radiogenic Neutron from the moderators:
```
/WLGD/detector/setGeometry baseline             # setting the geometry of the detector to the baseline design
/WLGD/event/saveAllEvents 0                     # only the Ge77 producing events are saved
/WLGD/detector/With_NeutronModerators 1         # using the moderator design with the tubes right around the re-entrance tubes
/WLGD/step/getDepositionInfo 1                  # save the information of multiplicity and total energy deposited in detectors
/run/initialize                                 
/WLGD/generator/setGenerator ModeratorNeutrons  # set the primary generator to the (Alpha,n) generator in the moderators
/run/beamOn 1000000
```

## Example for Ge77 production using Gd water and Turbine-like Moderators by Musun code:
```
/WLGD/detector/setGeometry baseline             # setting the geometry of the detector to the baseline design
/WLGD/event/saveAllEvents 0                     # only the Ge77 producing events are saved
/WLGD/detector/With_NeutronModerators 2         # using the moderator design with the tubes right around the re-entrance tubes
/WLGD/detector/With_Gd_Water 1                  # using the Gd in the water
/WLGD/detector/TurbineAndTube_Radius 200       # set the radius on which the center of mass of the pannels are alligned on [cm]
/WLGD/detector/TurbineAndTube_Length 100       # set the length of the moderator pannelss [cm]
/WLGD/step/getDepositionInfo 1                  # save the information of multiplicity and total energy deposited in detectors
/run/initialize                                 
/WLGD/generator/setGenerator Musun              # set the primary generator to the (Alpha,n) generator in the moderators
/WLGD/generator/setMUSUNFile path/to/file       # see the example/example_musun_file.dat
/run/beamOn 100                                 # should never exceed the size of the musun input file
```

## Example for investigating the output of all muons and their individual energy depositions in the whole cryostat
```
/WLGD/detector/setGeometry baseline             # setting the geometry of the detector to the baseline design
/WLGD/event/saveAllEvents 1                     # only the Ge77 producing events are saved
/WLGD/step/getDepositionInfo 1                  # save the information of multiplicity and total energy deposited in detectors
/WLGD/step/getIndividualDepositionInfo 1        # save the information of individual energy depositions inside the cryostat
/run/initialize                                 
/WLGD/generator/setGenerator Musun              # set the primary generator to the (Alpha,n) generator in the moderators
/WLGD/generator/setMUSUNFile path/to/file       # see the example/example_musun_file.dat
/run/beamOn 100                                 # should never exceed the size of the musun input file
```