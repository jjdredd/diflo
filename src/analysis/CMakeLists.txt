add_library(analysis OBJECT dist.cpp fft.cpp handedness.cpp momentum.cpp
  grid.cpp helicity.cpp)
target_include_directories(analysis PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(analysis -lm -lgsl -lgslcblas)
