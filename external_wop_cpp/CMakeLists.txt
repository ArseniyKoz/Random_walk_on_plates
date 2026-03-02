cmake_minimum_required(VERSION 3.20)

project(wop_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WOP_BUILD_TESTS "Build C++ tests" ON)

add_library(wop_core STATIC
    src/math/vec3.cpp
    src/geometry/polyhedron.cpp
    src/geometry/box.cpp
    src/rng/rng.cpp
    src/sampling/sampling.cpp
    src/estimation/estimation.cpp
    src/solver/wop_solver.cpp
    src/solver/wos_box_solver.cpp
)

target_include_directories(wop_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(wop_core PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(wop_core PRIVATE /W4 /permissive-)
else()
    target_compile_options(wop_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(wop_cli
    app/main.cpp
)
target_link_libraries(wop_cli PRIVATE wop_core)

if(WOP_BUILD_TESTS)
    enable_testing()

    add_executable(wop_core_smoke_tests
        tests/test_core_smoke.cpp
    )
    target_link_libraries(wop_core_smoke_tests PRIVATE wop_core)
    add_test(NAME wop_core_smoke_tests COMMAND wop_core_smoke_tests)

    add_executable(wop_sampling_stats_tests
        tests/test_sampling_stats.cpp
    )
    target_link_libraries(wop_sampling_stats_tests PRIVATE wop_core)
    add_test(NAME wop_sampling_stats_tests COMMAND wop_sampling_stats_tests)
endif()
