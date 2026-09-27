cmake_minimum_required(VERSION 3.24)
project(ember VERSION 0.1.0 LANGUAGES CXX)

option(EMBER_TWO_METAL_COLLISIONS
  "Approximate fixed GS98 collision carriers with two moment-matched representatives" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- Apple silicon -----------------------------------------------------------
# Target the native CPU.
# The small Henyey blocks use a pivoted kernel. -ffp-contract=fast
# permits FMA contraction in interpolation without the broader assumptions
# licensed by -ffast-math.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mcpu=apple-m4" EMBER_HAS_M4)
add_library(ember_flags INTERFACE)
target_compile_options(ember_flags INTERFACE
  $<$<CONFIG:Release>:-O3;-ffp-contract=fast;-fno-math-errno>
  $<$<CONFIG:Debug>:-O0;-g;-fno-omit-frame-pointer>
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
if(EMBER_HAS_M4)
  target_compile_options(ember_flags INTERFACE -mcpu=apple-m4)
elseif(APPLE)
  target_compile_options(ember_flags INTERFACE -mcpu=apple-m1)
endif()

# Deliberately NOT -ffast-math: it licenses the compiler to assume no NaN/Inf,
# and a stellar model legitimately probes states where a table lookup returns
# one.  We want those to surface, not to be optimised away.

find_package(Threads REQUIRED)
target_link_libraries(ember_flags INTERFACE Threads::Threads)

add_subdirectory(src)
add_subdirectory(apps)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
