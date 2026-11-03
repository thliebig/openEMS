# C++ unit tests, run with `ctest` in the build directory.
#
# Each test is a self-contained executable linking the openEMS library that was
# just built. Add a new one to TESTS below, no further wiring needed.

set(TESTS
  test_arraylib
  test_engine_plan
)

foreach(test ${TESTS})
  add_executable(${test} ${test}.cpp)
  # the tests use the headers of the source tree, not the installed ones
  target_include_directories(${test} PRIVATE ${openEMS_SOURCE_DIR})
  target_link_libraries(${test} openEMS)
  add_test(NAME ${test} COMMAND ${test})
endforeach()
