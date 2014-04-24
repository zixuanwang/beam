******************************************************************************
*   Beamformer
******************************************************************************
*   by Daniel Wang, 14th April 2014
******************************************************************************

******************************************************************************
*   Before compiling, make sure you have installed cmake and the compiler
*   supports c++11.
******************************************************************************

----------------------------------------------------------
Building the project using CMake from the command-line:
----------------------------------------------------------
Linux:
    mkdir build
    cd build
    cmake ../src
    make

Windows (MS Visual Studio):
    mkdir build
    cd build
    cmake ../src
    start beam.sln 
    
- A static library will be written to the "lib" directory.
- The execuables can be found in the "bin" directory.

----------------------------------------------------------
Running the various programs:
Help will be shown when you run each program without argument
----------------------------------------------------------
* On Linux or Mac:	./bin/program_name
* On Windows:		bin\Debug\program_name
