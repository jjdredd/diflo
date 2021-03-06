project(diflo)

# set CFLAGS etc

set(WFLAGS "-W -Wall -Wextra -std=c++11")
set(OPTFLAGS "-march=native -ffast-math -O3 -funroll-loops")
set(CMAKE_CXX_FLAGS  "${WFLAGS} ${OPTFLAGS}")

add_executable(diflo.elf diflo.cpp)
add_executable(sdnopt.elf sdnopt.cpp)
add_executable(multip.elf multip.cpp)

target_include_directories(diflo.elf PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(sdnopt.elf PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(multip.elf PUBLIC "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(input)
add_subdirectory(analysis)

# might actually need -lgslblas or something
target_link_libraries(diflo.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas)
target_link_libraries(sdnopt.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas)
target_link_libraries(multip.elf PUBLIC input
			     PUBLIC analysis -lm -lgsl -lgslcblas)
