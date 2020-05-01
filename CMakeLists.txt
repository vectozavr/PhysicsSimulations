cmake_minimum_required(VERSION 3.7)

project(PhysicsSimulations)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -static")
set(CMAKE_CXX_STANDARD 17)

add_library(vemath STATIC vemath.h vemath.cpp)

add_executable(TURTLE turtle.cpp settings.h)
add_executable(SALUTE salute.cpp settings.h)
add_executable(POTENTIAL potential.cpp settings.h)
add_executable(SEMICLASSICAL_OSCILLATION SemiclassicalOscillation.cpp gnuplot.h)
add_executable(LINEAR_DIFFERENTIAL_EQUATION linearDifferentialEquation.cpp gnuplot.h)
add_executable(DIFF_EQUATION_1 diffEquation1.cpp gnuplot.h)

add_executable(COVID19 covid19.cpp settings.h)

add_executable(DISCRETE_FOURIER discreteFouriertransform.cpp settings.h gnuplot.h)
add_executable(CONVOLUTION_AND_CROSS_CORELATION convolutionAndCrossCorelation.cpp settings.h gnuplot.h)
add_executable(LOW_FREQUENCY_FILTER lowFrequaencyFilter.cpp settings.h gnuplot.h)
add_executable(HIGH_FREQUENCY_FILTER highFrequencyFilter.cpp settings.h gnuplot.h)
add_executable(NOTCH_FILTER notchFilter.cpp settings.h gnuplot.h)
add_executable(SAMPLE_RATE_CONVERTION SampleRateConversion.cpp settings.h gnuplot.h)

# SFML Librares
set(SFML_DIR "D:/Libraries/SFML/lib/cmake/SFML")
set(SFML_STATIC_LIBRARIES TRUE)
find_package(SFML 2.5.1 COMPONENTS graphics audio REQUIRED)
include_directories(${SFML_INCLUDE_DIR})

target_link_libraries(TURTLE sfml-audio sfml-network sfml-graphics sfml-window sfml-system)
target_link_libraries(SALUTE sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)
target_link_libraries(POTENTIAL sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)
target_link_libraries(SEMICLASSICAL_OSCILLATION vemath)
target_link_libraries(LINEAR_DIFFERENTIAL_EQUATION vemath)
target_link_libraries(DIFF_EQUATION_1 vemath)

target_link_libraries(COVID19 sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)

target_link_libraries(DISCRETE_FOURIER vemath)
target_link_libraries(CONVOLUTION_AND_CROSS_CORELATION vemath)
target_link_libraries(LOW_FREQUENCY_FILTER vemath)
target_link_libraries(HIGH_FREQUENCY_FILTER vemath)
target_link_libraries(NOTCH_FILTER vemath)
target_link_libraries(SAMPLE_RATE_CONVERTION vemath)

# fftw
# target_link_libraries(main ${CMAKE_SOURCE_DIR}/libbingitup.a)

#target_link_libraries(CONVOLUTION_AND_CROSS_CORELATION ${CMAKE_SOURCE_DIR}/fftw/libfftw3-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3f-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3l-3.dll)
#target_link_libraries(HIGH_FREQUENCY_FILTER ${CMAKE_SOURCE_DIR}/fftw/libfftw3-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3f-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3l-3.dll)