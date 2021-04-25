# Detailed simulation of Forward Calorimeter Layer 1 
- Edited by Anson Kost, with the help of [Prof. John Rutherfoord](https://w3.physics.arizona.edu/people/john-rutherfoord)
- Author: Zhaoyuan "Maxwell" Cui, Advisor: [Prof. Erich Varnes](http://w3.physics.arizona.edu/people/erich-varnes)
- Physics Department, University of Arizona

## Running the Simulation

1. Build the program in a folder of your choice (see Maxwell's instructions below or the [Geant4 Book for Application Developers](https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/GettingStarted/makeFile.html)).
2. If you want to save data to files, then create a folder named "data".
3. In the parent folder of the data folder, enter the following command in the terminal (while in the parent of the data folder): `</path/to/FCal> </path/to/runBetas.mac>`. Here, `<path/to/FCal>` is the (relative) path to the FCal program in the build directory, and `</path/to/runBetas.mac>` is the (relative) path to the runBetas.mac script included in this repository.
4. Rename the "data" folder to save the data, and create a new, empty folder named "data" for the next run.

## Editing the Simulation

- To control the x and y offset of the position of the electron beam, edit the `offset` variable in src/FCalPrimaryGeneratorAction.cc.
- The energy of the electron and the spread of the beam can also be edited in src/FCalPrimaryGeneratorAction.cc.
- The source of the electrons can be switched between the external beam and the hot foil inside the electrodes by editing the file src/FCalActionalInitialization.cc.

Maxwell's original documentation is below:

# Detailed simulation of Forward Calorimeter Layer 1 ]
- Author: Zhaoyuan "Maxwell" Cui
- Advisor: Prof. [Erich Varnes](http://w3.physics.arizona.edu/people/erich-varnes)
- Physics Department, University of Arizona
- This is an undergraduate directed research

## Introduction
Forward Calorimeter (FCal) is one of the detectors in the ATLAS experiment located at Large Hadron Collider (LHC) CERN.
The first layer of FCal (FCal1) is a EM calorimeter.
Due to the large |η|, FCal1 is important in determining the location of the vertex of large pseudorapidity electrons 
and positrons.
The goal of the study is to use GEANT4 to simulate the FCal1 detector and the detecting process. By using the data from
the simulation, the timing resolution of FCal1 can be calculated. For a detail description of the simulation, please look
at the report [`Studying the timing resolution of FCal1 with a detailed simulation using Geant4.pdf`](https://github.com/maxwellcui/FCal1_GEANT4_simulation/blob/master/Studying%20the%20timing%20resolution%20of%20FCal1%20with%20a%20detailed%20simulation%20using%20Geant4.pdf)

## Before compiling
GEANT4 uses CMake to generate Makefile. Therefore, a CMakeLists and the propriate paths of files are important. 
The flags for cmake may vary on different machine and also under differnet circumstances of desired results. 

## Compile
In the directory of `CMakeLists.txt`, create a new directory called `build`:
```shell-script
mkdir build
```
Then change directory to the build directory:
```
cd build
```
Source Geant4 and invoke cmake:
```
source /PathToGeant4/bin/geant4.sh
cmake -DGeant4_DIR=/PathToGeant4/lib/Geant4-x.x.x ..
```
where the x's stand for the version of the Geant4.

A Makefile should be generated after the cmake process, to compile the program, type:
```
make 
```
or, if cpu has more than one core, one can use multiple cores to accelerate the compile process. Type:
```
make -jn
```
where n is the number of cores that the user want to use.

## Running program
After compiling the program, one can now run the program. Look for the executable in the `build` director. 
In the default `CMakeLists.txt`, macro files are copied to the `build` directory after successfully invoking cmake.
In the default setup of the program, running the program leads to the interactive mode of Geant4, which is a GUI that will
draw the geometry of the detector. In a detailed simulation, this can take significant amount of time and may lead to 
internal error if memory is not enough for the graphic drawing. 
For this reason, batch mode is recommended for collision simulation.
Intercative mode can be skipped by executing the program with macro file in an actual simulaton.
Nevertheless, in the phase of building the geometry of the detector, interactive mode is useful as it visiually shows the 
simulated geometry of the detector and can aid the user to create the desired geometry.

## Notes
Please note that Geant4 is highly customizable tool kit. User should design and/or modify the code to meet his or her own
concern.
Most parts, especially the action and the macro files, are therefore subject to change.
