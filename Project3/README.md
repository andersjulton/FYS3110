-----------main.cpp-----------
- Contains main program
- compile: make program
- run: ./main.exe 

-----------NSB.cpp-----------
- Contains vitrual class NBS for n-body systems
- Simulates position of all bodies in the System
	 using Euler's method and the velocity Verlet method

-----------SolarSystem.cpp-----------
- Contains class SolarSystem, child class of NSB 
- Contains class SolarSystemRelativistic, child class of SolarSystem 

-----------NSB.h-----------
- Header file for class NSB and children
- Contains struck MassObject

-----------planetsLib.h-----------
- Contains StellarObjectsLibrary

-----------utils.cpp-----------
- Contains general supporting functions