add_executable(ember-info info.cpp)
target_link_libraries(ember-info PRIVATE ember)
add_executable(ember-polytrope polytrope.cpp)
target_link_libraries(ember-polytrope PRIVATE ember)
