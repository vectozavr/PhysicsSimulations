cmake_minimum_required(VERSION 3.16)

project(PhysicsSimulations)
set(CMAKE_CXX_STANDARD 20)

add_library(vemath SHARED vemath.h vemath.cpp)

# semester 6

add_executable(TURTLE turtle.cpp settings.h)
add_executable(SALUTE salute.cpp settings.h)
add_executable(POTENTIAL potential.cpp settings.h)
add_executable(SEMICLASSICAL_OSCILLATION SemiclassicalOscillation.cpp gnuplot.h)
add_executable(LINEAR_DIFFERENTIAL_EQUATION linearDifferentialEquation.cpp gnuplot.h)
add_executable(DIFF_EQUATION_1 diffEquation1.cpp gnuplot.h)
add_executable(BROWNIAN_MOTION brownianMotion.cpp gnuplot.h)
add_executable(ANY_RANDOM_GRAPH AnyRandomGraph.cpp gnuplot.h)
add_executable(STATISTICS_IN_SILICON statisticsInSilicon.cpp gnuplot.h)

add_executable(POINT_FRAC pointsFraq.cpp)


add_executable(COVID19 covid19.cpp settings.h)

add_executable(DISCRETE_FOURIER discreteFouriertransform.cpp settings.h gnuplot.h)
add_executable(CONVOLUTION_AND_CROSS_CORELATION convolutionAndCrossCorelation.cpp settings.h gnuplot.h)
add_executable(LOW_FREQUENCY_FILTER lowFrequaencyFilter.cpp settings.h gnuplot.h)
add_executable(HIGH_FREQUENCY_FILTER highFrequencyFilter.cpp settings.h gnuplot.h)
add_executable(NOTCH_FILTER notchFilter.cpp settings.h gnuplot.h)
add_executable(SAMPLE_RATE_CONVERTION SampleRateConversion.cpp settings.h gnuplot.h)

# semester 7

add_executable(MACHINE_EPSILON machine_epsilon.cpp)
add_executable(SCHRODINGER_EQUATION schrodinger_equation.cpp)
add_executable(CALCULATE_INTEGRALS calculate_integrals.cpp)

add_executable(BESSELS_FUNCTION besselsFunction.cpp)
add_executable(INTROPOLATION_POLINOM interpolationPolinom.cpp)
add_executable(KOSHI_PROBLEM koshiProblem.cpp)
add_executable(PREDATOR_PREY predatorPrey.cpp)

add_executable(RIGID_SYSTEM rigidSystem.cpp)
add_executable(SWEEP_METHOD sweepMethod.cpp)

#Crank–Nicolson method Inverse iteration Spectral density
add_executable(CRANK_NICOLSON_METHOD crankNicolsonMethod.cpp)
add_executable(INVERSE_ITERATION inverseIteration.cpp)
add_executable(SPECTRAL_DENSITY spectralDensity.cpp)

add_executable(COURSE_TASK courseTask.cpp)

# SFML Librares
##
if(${APPLE})
    include_directories(/usr/local/include)
else()

    set(SFML_DIR "D:/Libraries/SFML/lib/cmake/SFML")
    set(SFML_STATIC_LIBRARIES TRUE)

endif()

find_package(SFML 2.5.1 COMPONENTS graphics audio REQUIRED)

if (SFML_FOUND)
    include_directories(${SFML_INCLUDE_DIR})
endif()
##

target_link_libraries(TURTLE sfml-audio sfml-network sfml-graphics sfml-window sfml-system)
target_link_libraries(SALUTE sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)
target_link_libraries(POTENTIAL sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)
target_link_libraries(SEMICLASSICAL_OSCILLATION vemath)
target_link_libraries(LINEAR_DIFFERENTIAL_EQUATION vemath)
target_link_libraries(DIFF_EQUATION_1 vemath)

target_link_libraries(COVID19 sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)
target_link_libraries(BROWNIAN_MOTION sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)

target_link_libraries(POINT_FRAC sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)



target_link_libraries(DISCRETE_FOURIER vemath)
target_link_libraries(CONVOLUTION_AND_CROSS_CORELATION vemath)
target_link_libraries(LOW_FREQUENCY_FILTER vemath)
target_link_libraries(HIGH_FREQUENCY_FILTER vemath)
target_link_libraries(NOTCH_FILTER vemath)
target_link_libraries(SAMPLE_RATE_CONVERTION vemath)
target_link_libraries(ANY_RANDOM_GRAPH vemath)


target_link_libraries(CALCULATE_INTEGRALS vemath)

target_link_libraries(BESSELS_FUNCTION vemath)
target_link_libraries(INTROPOLATION_POLINOM vemath)
target_link_libraries(KOSHI_PROBLEM vemath)

target_link_libraries(PREDATOR_PREY vemath)
target_link_libraries(RIGID_SYSTEM vemath)
target_link_libraries(SWEEP_METHOD vemath)

target_link_libraries(CRANK_NICOLSON_METHOD sfml-audio sfml-network sfml-graphics sfml-window sfml-system vemath)
target_link_libraries(INVERSE_ITERATION vemath)
target_link_libraries(SPECTRAL_DENSITY vemath)
target_link_libraries(COURSE_TASK vemath)

target_link_libraries(STATISTICS_IN_SILICON vemath)

# fftw
# target_link_libraries(main ${CMAKE_SOURCE_DIR}/libbingitup.a)

#target_link_libraries(CONVOLUTION_AND_CROSS_CORELATION ${CMAKE_SOURCE_DIR}/fftw/libfftw3-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3f-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3l-3.dll)
#target_link_libraries(HIGH_FREQUENCY_FILTER ${CMAKE_SOURCE_DIR}/fftw/libfftw3-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3f-3.dll ${CMAKE_SOURCE_DIR}/fftw/libfftw3l-3.dll)