# warwick-legend
Legend muon background simulations

This is a LEGEND Monte-Carlo simulation application,
suitable for large production runs estimating Germanium-77 production in
Germanium-76 crystals in different underground laboratories.

The ItalSim branch of warwick-legend contains updates to expand the functionality
of the simulation. In addition to muon studies, it is also possible to perform
dedicated neutron studies, nuclear decay studies, and optical photon
scintillation and detection studies.


## Requirements

Multi-threading operation and MPI capability

ROOT ntuple output

Operate with command line input and Geant4 macros.

Macro commands control many aspects of the simulation, including primary
particles, geometry options, and optical implementation.

CLI controls number of threads, macro file, output file name (for production runs)

User limits can be set for runtime optimization


## How to Build/Develop
The project has the following requirements:

- Linux/macOS system (only CentOS7, Catalina tested at present)
- C++17 compatible compiler (GCC9, Xcode 11 tested at present)
- CMake 3.12 or newer
- Geant4 10.7 built with multithreading and gdml support
  - Qt or OpenGL/X11 driver support if you want to visualize the geometry/tracks/hits

It may be compiled using:

```console
$ mkdir build && cd build
$ cmake ..
$ make
```

The resulting `warwick-legend` application may be run without arguments to start an interactive
session. Otherwise run `warwick-legend --help` to see a list of options for batch mode running.

In order to run the simulation, call the executable which was compiled with 'make'

```
$/path/to/build/warwick-legend
```

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
for fixes will be output during building. At present, the `-fix` option to
automatically apply the suggested changes is not used in order to leave the
decision up to the developer. The set of fixes applied are:

- `readability-*`
- `modernize-*`
- `performance-*`

For a full listing of the wildcards, see [the `clang-tidy` documentation](https://clang.llvm.org/extra/clang-tidy/checks/list.html).

## Code Details
### Cross section bias reference
M.H. Mendenhall and R.A. Weller, NIM A667 (2012) 38-43

### Ntuple output columns (in the "Score" tree)
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

A lot has been added in my version of the fork. - Moritz
An output for the Ge77m veto has been added, more options to adjust the geometry via macros and several new primary generators.

## Overview over the Macros

### Runaction Macro
Macros regarding the output of the simulation
All output macros have 0 as the default and must be set to 1 for the relevant output to fill
```

Macros added by Moritz:
/WLGD/runaction/
  - WriteOutNeutronProductionInfo
  - WriteOutGeneralNeutronInfo
  - WriteOutAllNeutronInfoRoot
  - WriteOutAdvancedMultiplicity
  - getIndividualGeDepositionInfo
  - getIndividualGdDepositionInfo
  - getNeutronCaptureSiblings
  - readMuonCrossingWLSR

Runaction macros added by CJ:
/WLGD/runaction/
  - WriteOpticalProductionData (beginning-of-track info about all optical photons produced, saved in "OpticksTracks" tree)
  - WriteOpticalMapData        (for generating the optical map, saved in "OpticalMapData" tree; probably only used by CJ)
  - WriteStepData              (for writing out all step-level info, saved in "Steps" tree; user adds their own analysis cuts in the SteppingAction.cc)
```
### Event Macro
Macro to adjust the condition to save all events (1) or just the ones with Ge77 production (0)
Only works with the "Score" output tree.
```
/WLGD/event/
  - saveAllEvents (no: [0], yes: 1)
```

Similarly, it is possible to reduce the Steps tree to fewer parameters with the command
```
/WLGD/runaction/ReduceStepData 1
```

### Generator Macro
Macros to control the primary generator 
```
/WLGD/generator/
  - depth
  - setMUSUNFile (path to individual MUSUN file)
  - setMUSUNDirectory (full path to directory containing MUSUN files)
  - setGenerator (options: "MeiAndHume", "Musun", "Ge77m", "Ge77andGe77m", "ModeratorNeutrons", "ExternalNeutrons", "SimpleGammaGun", "SimpleNeutronGun", "OpticalMap", "ArgonCaptureGammas")
```

More information on SetMUSUNDirectory can be found in the OpenMUSUNDirectory method of src/WLGDPrimaryGeneratorAction.cc

### Detector Macro
Macros to control the detector geometry
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
- PolygonShield (for flat-sided neutron shield only, option 5 in With_NeutronModerators)
    - NSides

###Optics commands (only the first is currently implemented)

/WLGD/optics/
- WithOptics          (1 for on, 0 for off)
- CladdingLayers      (default is 1)
- CladdingMaterial    (default PEN)
- CladdingThickness   (in um, default is 100)
- LightGuideLength    (in cm, default 100)
- LightGuideWidth     (default 3 cm)
- LightGuideMaterial  (default PMMA)
- LightGuideSpacing   (default is 30 cm)
- NLightGuides        (default is 12)
- UseWLSCoating       (1 for on, 0 for off)
- WLSCoatingMaterial  (default is TPB)

If XeConc > 0, the LAr optical properties will automatically change to the properties of LAr doped with 100 ppm xenon. Currently, this is the only xenon doping concentration implemented.

For the number of light guides and light guide spacing, the user can choose to set both. However, if only one is provided, a reasonable value for the other will be calculated based on the shield height.


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
Macros to adjust whether additional output (additional to the Ge77 production) is recorded in the first place ("Score" tree only)
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
/WLGD/generator/setMUSUNDirectory path/to/directory #Use a directory containing any number of MUSUN files, chosen randomly
/run/beamOn 100                                 # 
```