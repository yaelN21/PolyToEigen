cmake_minimum_required(VERSION 3.2.2)
project(Triangulation)
set (CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)
find_package (Eigen3 REQUIRED)
message(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include/)
enable_testing()
add_library(Triangulation ${CMAKE_CURRENT_SOURCE_DIR}/src/Poly.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/PolyBase.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/TriangulationBase.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/LinearLS.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/TriangulationBase.cpp)
target_link_libraries(Triangulation ${OpenCV_LIBS} Eigen3::Eigen )
add_executable(main main.cpp)
target_link_libraries(main Triangulation)


#add_test(AllTests Triangulation)

#configure_file(TriangulationConfig.cmake ${CMAKE_INSTALL_PREFIX}/lib/cmake/Triangulation/TriangulationConfig.cmake)

