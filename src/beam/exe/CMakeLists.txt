ADD_EXECUTABLE(test_ssl test_ssl.cpp)

TARGET_LINK_LIBRARIES(test_ssl libbeam ${CMAKE_PREFIX_PATH}/libfftw3f-3.lib ${CMAKE_PREFIX_PATH}/libfftw3-3.lib ${CMAKE_PREFIX_PATH}/libfftw3l-3.lib)