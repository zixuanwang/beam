******************************************************************************
*   Geometric VLAD
******************************************************************************
*   by Zixuan Wang, 14th September 2013
*   http://www.stanford.edu/~zxwang
******************************************************************************

******************************************************************************
*   Before compiling, make sure you have installed boost, cmake and opencv
******************************************************************************

----------------------------------------------------------
Building the project using CMake from the command-line:
----------------------------------------------------------
Linux:
    export OpenCV_DIR="~/OpenCV/build"
    export Boost_ROOT="~/boost"
    mkdir build
    cd build
    cmake -D OpenCV_DIR=$OpenCV_DIR -D Boost_ROOT=$Boost_ROOT ../src
    make

Windows (MS Visual Studio):
    set OpenCV_DIR="C:\opencv\build"
    set Boost_ROOT="C:\boost"
    mkdir build
    cd build
    cmake -G "Visual Studio 11 Win64" -D OpenCV_DIR=%OpenCV_DIR% -D Boost_ROOT=%Boost_ROOT% ../src
    start gvlad.sln 
    
- A static library will be written to the "lib" directory.
- The execuables can be found in the "bin" directory.

----------------------------------------------------------
Running the various programs:
Help will be shown when you run each program without argument
----------------------------------------------------------
* On Linux or Mac:	./bin/program_name
* On Windows:		bin\Debug\program_name

----------------------------------------------------------
Example to use provided executables:
----------------------------------------------------------


