add_library(ember
  composition.cpp
  fermi.cpp
  interp.cpp
  opacity_ferguson.cpp
  eos_components.cpp
  eos_ideal.cpp)
target_include_directories(ember PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(ember PUBLIC ember_flags)
