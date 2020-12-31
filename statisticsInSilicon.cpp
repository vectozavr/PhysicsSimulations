//
// Created by Иван Ильин on 28.11.2020.
//

#include <cmath>
#include <iostream>
#include "gnuplot.h"
#include "vemath.h"

using namespace std;
using namespace vemath;

class Silicon {
    private:
        /*
        // Физические константы:
        double h = 1.0f;
        double k = 1.0f;
        */
        double h = 6.62607015e-34;
        double k = 1.380649e-23;

        /*
         * Параметры для ввода в модель ( http://www.ioffe.ru/SVA/NSM/Semicond/Si/electric.html#Basic ):
         * 1) E_d - положение уровня донора
         * 2) E_g - ширина запрещённой зоны
         * 3) E_c - дно зоны проводимости
         * 4) m - эффективная масса носителей заряда
         * 5) N_d0 - концентрация доноров (задаётся от 10^15 до 10^22 на cm^3)
         * 6) T_0 - начальная температура
         * 7) T_1 - конечная температура
         * Все энергии отсчитывает от потолка валентной зоны, только энергия донора вниз от дна зоны проводимости.
         * TODO: Программа переводит все в единицы СГС (или в СИ по желанию).
        */

        double E_d;
        double E_g;
        double E_c;
        double m;
        double N_d0;

        ComplexPlot F_T; // dependence F on T
        ComplexPlot n_T; // dependence n on T

        ComplexPlot eq_F; // dependence eq on F

        double f_mult_g(double E, double F, double T) {
            return 1.0f/(2.0f*PI*PI)*pow(2.0f*m/(h*h), 1.5f) * sqrt(E - E_c)/(exp((E - F) / (k*T)) + 1);
        }

        double eq(double F, double T, int Nn = 1000) {
            return E_g - E_d - k*T*log(N_d0/n(F, T, Nn) - 1) - F;
        }

        double deq(double F, double T, double dF = 0.0001f, int Nn = 1000){
            double eq2 = eq(F + dF, T, Nn);
            double eq1 = eq(F - dF, T, Nn);

            return (eq2 - eq1) / (2.0f * dF);
        }

    public:
        Silicon(double E_d, double E_g, double E_c, double m, double N_d0)
        {
            this->E_d = E_d;
            this->E_g = E_g;
            this->E_c = E_c;
            this->m = m;
            this->N_d0 = N_d0;
        }

        double n(double F, double T, int N = 1000) {
            double inf = 1000.0f * E_c;

            double dE = (inf - E_c) / N;
            double E = E_c;

            double S = 0.0f;

            while (E < inf) {
               S += f_mult_g(E, F, T) + 4.0f*f_mult_g(E + dE/2.0f, F, T) + f_mult_g(E + dE, F, T);
               E += dE;
            }

            return S * dE / 6.0f;
        }

        void calcilate_F_from_T(double T_0, double T_1, double tol = 0.001f, int NT = 1000, int Nn = 1000){
            F_T.clear();
            n_T.clear();

            double dT = (T_1 - T_0) / NT;

            double T = T_0;
            while(T < T_1) {
                // Метод Ньютона
                double F_N = 0.0f; // искомое F
                double F1_N = 0.7f; // нулевое приближение
                double dF = 0.0001f; // шаг дифференцирования
                int s = 0; // число итераций

                while(abs(F_N - F1_N) > tol) {
                    F_N = F1_N;
                    F1_N = F_N - eq(F_N + dF, T, Nn) / deq(F_N + dF, T, dF, Nn);
                    s++;
                }
                F_T.push(T, F_N);
                n_T.push(T, n(F_N, T, Nn));

                T += dT;
            }
        }

        void plot_eq_from_F(double F0, double F1, double T, int NF = 1000, int Nn = 1000) {
            eq_F.clear();
            double F = F0;
            double dF = (F1 - F0) / NF;

            while(F < F1) {
                eq_F.push(F, eq(F, T, Nn));
                F += dF;
            }

            saveVectorPoint2DToFile(eq_F.real(), "eq_F.dat");
            GnuplotPipe gp;
            gp.sendLine(R"(plot "eq_F.dat" with lines)");
        }

        bool saveData() {
            if(saveVectorPoint2DToFile(F_T.real(), "F_T.dat") && saveVectorPoint2DToFile(n_T.real(), "n_T.dat"))
                return true;
            else
                return false;
        }

        void plotData() {
            GnuplotPipe gp;
            gp.sendLine(R"(set multiplot layout 2,1)");
            gp.sendLine(R"(plot "F_T.dat" with lines)");
            gp.sendLine(R"(plot "n_T.dat" with lines)");
        }

        void plotPNGData() {
            GnuplotPipe gp;
            gp.sendLine(R"(set term png)");
            gp.sendLine(R"(set output "f_from_t.png")");
            gp.sendLine(R"(plot "F_T.dat" with lines)");
            gp.sendLine(R"(unset output)");

            gp.sendLine(R"(set term png)");
            gp.sendLine(R"(set output "n_from_t.png")");
            gp.sendLine(R"(plot "n_T.dat" with lines)");
            gp.sendLine(R"(unset output)");
        }

};

int main() {

    /*
     * Параметры для ввода в модель ( http://www.ioffe.ru/SVA/NSM/Semicond/Si/electric.html#Basic ):
     * 1) E_d - положение уровня донора
     * 2) E_g - ширина запрещённой зоны
     * 3) E_c - дно зоны проводимости
     * 4) m - эффективная масса носителей заряда
     * 5) N_d0 - концентрация доноров (задаётся от 10^15 до 10^22 на cm^3)
     * 6) T_0 - начальная температура
     * 7) T_1 - конечная температура
     * Все энергии отсчитывает от потолка валентной зоны, только энергия донора вниз от дна зоны проводимости.
     * TODO: Программа переводит все в единицы СГС (или в СИ по желанию).
    */

    /*
    double E_d = 1.0f;
    double E_g = 1.0f;
    double E_c = 1.0f;
    double m = 1.0f;
    double N_d0 = 1.0f;
    double T_0 = 0.0f;
    double T_1 = 100.0f;
    */
    double E_d = 2.0f * 1.60e-19;
    double E_g = 1.5f * 1.60e-19;
    double E_c = 1.5f * 1.60e-19;
    double m = 9.1093837015e-31;
    double N_d0 = 1e22;
    double T_0 = 0.0f;
    double T_1 = 100.0f;

    // TODO: solve main equation: F = E_g - E_d - k*T * ln(N_d0/n(F) - 1)
    int Nn = 1000;
    int NT = 1000;
    double tol = 0.001f;

    Silicon silicon(E_d, E_g, E_c, m, N_d0);
    silicon.calcilate_F_from_T(T_0, T_1, tol, NT, Nn);
    silicon.saveData();

    silicon.plotPNGData();
    silicon.plotData();

    //silicon.plot_eq_from_F(-60, 60, 30);
}
