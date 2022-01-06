cmake_minimum_required(VERSION 3.18)

project(diflo)

find_package(Boost 1.36.0 REQUIRED COMPONENTS program_options)

# set CFLAGS etc
set(WFLAGS "-W -Wall -Wextra -g -std=c++11")
set(OPTFLAGS "-march=native -ffast-math -O3 -funroll-loops")
set(CMAKE_CXX_FLAGS  "${WFLAGS} ${OPTFLAGS}")

add_executable(diflo.elf diflo.cpp)
add_executable(sdnopt.elf sdnopt.cpp)
add_executable(multip.elf multip.cpp)
add_executable(sphand.elf sphand.cpp)

include_directories(AFTER ${Boost_INCLUDE_DIR})
include_directories(AFTER "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(input)
add_subdirectory(analysis)
add_subdirectory(test)

# might actually need -lgslblas or something
target_link_libraries(diflo.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas)
target_link_libraries(sdnopt.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas)
target_link_libraries(multip.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas)

target_link_libraries(sphand.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas PUBLIC ${Boost_LIBRARIES})
