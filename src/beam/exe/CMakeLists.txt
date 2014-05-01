ADD_EXECUTABLE(beamformer beamformer.cpp)
ADD_EXECUTABLE(test test.cpp)

TARGET_LINK_LIBRARIES(beamformer libbeam)
TARGET_LINK_LIBRARIES(test libbeam)