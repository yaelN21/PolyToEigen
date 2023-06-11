cmake_minimum_required(VERSION 3.2.2)
project(Triangulation)
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -pedantic -Wextra -Wconversion -march=native -Wall -I. -ftest-coverage -fprofile-arcs")
set(Eigen3_DIR "/usr/share/eigen3/cmake")
find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
enable_testing()


set (HEADERS "Poly.h" "PolyBase.h" "TriangulationBase.h" "LinearLS.h")
set (SOURCES "Poly.cpp" "PolyBase.cpp" "TriangulationBase.cpp" "LinearLS.cpp")

add_executable(Triangulation main.cpp ${SOURCES})

target_link_libraries(Triangulation  Eigen3::Eigen)


add_test(AllTests Triangulation)

add_library(triangulation SHARED ${SOURCES})
set_target_properties(triangulation PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(triangulation PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")

install(TARGETS triangulation LIBRARY DESTINATION lib/Triangulation INCLUDES DESTINATION include/Triangulation)
install(FILES ${HEADERS} DESTINATION include/Triangulation)
configure_file(TriangulationConfig.cmake ${CMAKE_INSTALL_PREFIX}/lib/cmake/Triangulation/TriangulationConfig.cmake)

