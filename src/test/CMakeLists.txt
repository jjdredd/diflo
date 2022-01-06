include_directories(AFTER ${Boost_INCLUDE_DIR})
include_directories(AFTER "${CMAKE_SOURCE_DIR}/include")

add_executable(rp_test.elf rp_test.cpp)
target_link_libraries(rp_test.elf PUBLIC input PUBLIC analysis -lm -lgsl -lgslcblas)

add_executable(stats_test.elf stats_test.cpp)
target_link_libraries(stats_test.elf PUBLIC input PUBLIC analysis -lm -lgsl -lgslcblas)
