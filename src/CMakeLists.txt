add_library(ember
  composition.cpp
  fermi.cpp
  interp.cpp
  opacity_table.cpp
  opacity_blend.cpp
  eos_components.cpp
  eos_cms19.cpp
  eos_helmholtz.cpp
  eos_composition.cpp
  nuclear_pp.cpp
  convection.cpp
  atmosphere_grey.cpp
  atmosphere_table.cpp
  boundary.cpp
  henyey.cpp
  relaxation.cpp
  evolution.cpp
  structure.cpp
  eos_ideal.cpp)
target_include_directories(ember PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(ember PUBLIC ember_flags)
