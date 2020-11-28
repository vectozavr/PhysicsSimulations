//
// Created by Иван Ильин on 28.11.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

double n(double F, double T, double E_c, double m, double h, double k, int N = 1000) {

    double inf = 100.0f * E_c;

    double dE = (inf - E_c) / N;
    double E = E_c;

    double S = 0.0f;

    while (E < inf) {
        //S += f(x) + 4.0f*f(x + h/2.0f) + f(x + h);
        S +=    sqrt(E - E_c)/(exp((E - F) / k*T) + 1) +
                4.0f*sqrt(E + dE/2.0f - E_c)/(exp((E + dE/2.0f - F) / k*T) + 1) +
                sqrt(E + dE - E_c)/(exp((E + dE - F) / k*T) + 1);
        E += dE;
    }

    return 1.0f/(2.0f*PI*PI)*pow(2.0f*m/(h*h), 1.5f) * S * dE / 6.0f;
}

int main() {

    // Физические константы:
    double h = 1.0f;
    double k = 1.0f;

    /*
     * Параметры для ввода в модель ( http://www.ioffe.ru/SVA/NSM/Semicond/Si/electric.html#Basic ):
     * 1) E_d - положение уровня донора
     * 2) E_g - запрещённая зона
     * 3) E_c - дно зоны проводимости
     * 3) m - эффективная масса носителей заряда
     * 4) N_d0 - концентрация доноров (задаётся от 10^15 до 10^22 на cm^3)
     * 5) T_0 - начальная температура
     * 6) T_1 - конечная температура
     * Все энергии отсчитывает от потолка валентной зоны, только энергия донора вниз от дна зоны проводимости.
     * TODO: Программа переводит все в единицы СГС (или в СИ по желанию).
    */

    double E_d = 1.0f;
    double E_g = 1.0f;
    double m = 1.0f;
    int N_d0 = 10;
    double T_0 = 10.0f;
    double T_1 = 100.0f;

    ComplexPlot F_T; // dependence F on T
    ComplexPlot n_T; // dependence n on T

    // TODO: solve main equation: F = E_g - E_d - k*T * ln(N_d0/n(F) - 1)


    saveVectorPoint2DToFile(F_T.real(), "F_T.dat");
    saveVectorPoint2DToFile(n_T.real(), "n_T.dat");

    GnuplotPipe gp;
    gp.sendLine(R"(set multiplot layout 2,1)");
    gp.sendLine(R"(plot "F_T.dat" with lines)");
    gp.sendLine(R"(plot "n_T.dat" with lines)");
}
