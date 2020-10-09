//
// Created by Иван Ильин on 17.09.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;

int main() {

    float   epsilon_float  = 3.0f;
    double  epsilon_double = 3.0f;

    while (1.0f + epsilon_float / 2.0f != 1.0f) {
        epsilon_float /= 2.0f;
    }

    while (1.0f + epsilon_double / 2.0f != 1.0f) {
        epsilon_double /= 2.0f;
    }

    cout << "--------------------------------------------" << endl;

    cout << "Epsilon float  = " << epsilon_float << endl;
    cout << "Epsilon double = " << epsilon_double << endl;

    cout << "--------------------------------------------" << endl;

    cout << "1 : 1 + E/2 : 1 + E : 1 + E + E/2 (for float)" << endl;
    printf("%20.16f ", 1.0f);
    printf("%20.16f ", 1.0f + epsilon_float / 2.0f);
    printf("%20.16f ", 1.0f + epsilon_float);
    printf("%20.16f ", 1.0f + epsilon_float + epsilon_float / 2.0f);

    cout << endl << "1 : 1 + E/2 : 1 + E : 1 + E + E/2 (for double)" << endl;
    printf("%20.20f ", 1.0f);
    printf("%20.20f ", 1.0f + epsilon_double / 2.0f);
    printf("%20.20f ", 1.0f + epsilon_double);
    printf("%20.20f ", 1.0f + epsilon_double / 2.0f + epsilon_double);
}