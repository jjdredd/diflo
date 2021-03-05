add_library(input OBJECT data.cpp)
target_include_directories(input PUBLIC "${CMAKE_SOURCE_DIR}/include")
