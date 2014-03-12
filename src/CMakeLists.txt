PROJECT(beam)

cmake_minimum_required(VERSION 2.8)

INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
)
SET( LIBRARY_OUTPUT_PATH ${beam_BINARY_DIR}/../lib CACHE PATH
     "Output directory for the gvlad library" )
SET( EXECUTABLE_OUTPUT_PATH 
     ${beam_BINARY_DIR}/../bin CACHE PATH
     "Output directory for the gvlad executables" )

SUBDIRS(beam/lib)
SUBDIRS(beam/exe)
