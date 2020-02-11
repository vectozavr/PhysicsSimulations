cmake_minimum_required(VERSION 3.7)

project(PhysicsSimulations)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -static")
set(CMAKE_CXX_STANDARD 17)


add_executable(TURTLE turtle.cpp settings.h)
add_executable(DISCRETE_FOURIER discreteFouriertransform.cpp settings.h gnuplot.h)

# SFML Librares
set(SFML_DIR "D:/Libraries/SFML/lib/cmake/SFML")
set(SFML_STATIC_LIBRARIES TRUE)
find_package(SFML 2.5.1 COMPONENTS graphics audio REQUIRED)
include_directories(${SFML_INCLUDE_DIR})

target_link_libraries(TURTLE sfml-audio sfml-network sfml-graphics sfml-window sfml-system)