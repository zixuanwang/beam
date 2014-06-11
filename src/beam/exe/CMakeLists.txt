ADD_EXECUTABLE(beamformer beamformer.cpp)

TARGET_LINK_LIBRARIES(beamformer libbeam armadillo blas lapack)
