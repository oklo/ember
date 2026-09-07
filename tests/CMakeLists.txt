add_executable(test_eos test_eos.cpp)
target_link_libraries(test_eos PRIVATE ember)
add_test(NAME eos COMMAND test_eos)
