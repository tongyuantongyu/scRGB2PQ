cmake_minimum_required(VERSION 3.20)
project(num_exp C CXX)

set(CMAKE_CXX_STANDARD 20)

set(SC2PQ_LEAN_MODE OFF CACHE BOOL "Use Lean implementation (Only generate code suitable for upto 12 bit precision)?")
set(SC2PQ_LOCAL_HWY ON CACHE BOOL "Use local highway library?")
set(SC2PQ_ENABLE_TEST OFF CACHE BOOL "Enable test?")
set(SC2PQ_LOCAL_GTEST ON CACHE BOOL "Use local googletest library?")

if (SC2PQ_LOCAL_HWY)
    set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "")
    set(HWY_ENABLE_INSTALL OFF CACHE BOOL "")
    set(HWY_ENABLE_TESTS OFF CACHE BOOL "")
    set(HWY_FORCE_STATIC_LIBS ON CACHE BOOL "")
    add_subdirectory(highway EXCLUDE_FROM_ALL)
else()
    find_package(libhwy REQUIRED)
endif()

add_library(sc2pq sc2pq.h sc2pq.cpp)
target_include_directories(sc2pq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sc2pq PUBLIC hwy)

# This special target is for checking compiler-generated assembly.
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC" OR CMAKE_CXX_SIMULATE_ID MATCHES "MSVC")
    add_library(sc2pq_asm OBJECT sc2pq_asm.cpp)
    set_source_files_properties(sc2pq_asm.cpp PROPERTIES
            COMPILE_OPTIONS "/FAu;/Fa${CMAKE_CURRENT_BINARY_DIR}\\sc2pq.asm")
    target_include_directories(sc2pq_asm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(sc2pq_asm PRIVATE hwy)
    set_target_properties(sc2pq_asm PROPERTIES EXCLUDE_FROM_ALL ON)
else ()
    add_library(_sc2pq_asm OBJECT sc2pq_asm.cpp)
    set_source_files_properties(sc2pq_asm.cpp PROPERTIES COMPILE_OPTIONS -S)
    target_include_directories(_sc2pq_asm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(_sc2pq_asm PRIVATE hwy)
    set_target_properties(_sc2pq_asm PROPERTIES EXCLUDE_FROM_ALL ON)
    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/sc2pq.asm
            COMMAND cmake -E copy $<TARGET_OBJECTS:_sc2pq_asm> ${CMAKE_CURRENT_BINARY_DIR}/sc2pq.asm
            DEPENDS _sc2pq_asm
            VERBATIM)
    add_custom_target(sc2pq_asm DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/sc2pq.asm)
endif ()



if (SC2PQ_LEAN_MODE)
    target_compile_definitions(sc2pq PRIVATE SC2PQ_LEAN)
    target_compile_definitions(sc2pq_asm PRIVATE SC2PQ_LEAN)
endif()

if (SC2PQ_ENABLE_TEST)
    enable_testing()

    if (SC2PQ_LOCAL_HWY)
        set(BUILD_GMOCK OFF CACHE BOOL "")
        set(INSTALL_GTEST OFF CACHE BOOL "")
        add_subdirectory(googletest EXCLUDE_FROM_ALL)
    else()
        find_package(GTest REQUIRED)
        add_library(gtest_main ALIAS GTest::gtest_main)
    endif()

    add_executable(sc2pq_test sc2pq_test.cpp)
    target_link_libraries(sc2pq_test PUBLIC sc2pq gtest_main)
endif()
