AUX_SOURCE_DIRECTORY("${PROJECT_SOURCE_DIR}/test" DIR_TESTS)
SET(SRC "${PROJECT_SOURCE_DIR}/src" )

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin/test")

add_executable(PolyInterp_test
  ${SRC}/PolyInterp.cpp
  ${DIR_TESTS}
)

target_link_libraries(
  PolyInterp_test
  cppunit
  ${CMAKE_DL_LIBS}
)

ADD_TEST(NAME PolyInterp_test COMMAND "${PROJECT_SOURCE_DIR}/bin/test/PolyInterp_test")
set_tests_properties (PolyInterp_test
PROPERTIES FAIL_REGULAR_EXPRESSION "failed")
