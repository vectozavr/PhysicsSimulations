//
// Created by ivan- on 11.02.2020.
//
#include <algorithm>
#include "gnuplot.h"
#include "vemath.h"

using namespace vemath;
using namespace std;

typedef complex<double> w_type;

static vector<w_type> fft(const vector<w_type> &In)
{
    int i = 0, wi = 0;
    int n = In.size();
    vector<w_type> A(n / 2), B(n / 2), Out(n);
    if (n == 1) {
        return vector<w_type>(1, In[0]);
    }
    i = 0;
    copy_if( In.begin(), In.end(), A.begin(), [&i] (w_type e) {
        return !(i++ % 2);
    } );
    copy_if( In.begin(), In.end(), B.begin(), [&i] (w_type e) {
        return (i++ % 2);
    } );

    vector<w_type> At = fft(A);
    vector<w_type> Bt = fft(B);

    transform(At.begin(), At.end(), Bt.begin(), Out.begin(), [&wi, &n]
            (w_type& a, w_type& b) {
        return  a + b * exp(w_type(0, 2 * M_PI * wi++ / n));
    });
    transform(At.begin(), At.end(), Bt.begin(), Out.begin() + n / 2, [&wi, &n]
            (w_type& a, w_type& b) {
        return  a + b * exp(w_type(0, 2 * M_PI * wi++ / n));
    });
    return Out;
}

void dft1(const ComplexPlot& data, ComplexPlot& transform)
{
    double pi2 = 2.0 * M_PI;
    double angleTerm,cosineA,sineA;
    double invs = 1.0 / data.size();
    for(unsigned int y = 0; y < data.size(); y++) {
        std::complex<double> res = {0, 0};
        for(unsigned int x = 0; x < data.size(); x++) {
            angleTerm = pi2 * y * x * invs;
            cosineA = cos(angleTerm);
            sineA = sin(angleTerm);
            res += std::complex<double>{data.v_c[x].second.real() * cosineA - data.v_c[x].second.imag() * sineA, 0};
            res += std::complex<double>{0, data.v_c[x].second.real() * sineA + data.v_c[x].second.imag() * cosineA};
        }
        res *= invs;
        transform.push(y, res);
    }
}
/*
    int ln = (int)floor( log(argc - 1.0) / log(2.0) );
    vector<w_type> In(1 << ln);
    std::transform(argv + 1, argv + argc, In.begin(),[&](const char* arg) {
        return w_type(atof(arg), 0);
    });
    vector<w_type> Out = fft(In);
    for (vector<w_type>::iterator itr = Out.begin(); itr != Out.end(); itr++) {
        cout  << *itr << endl;
    }
*/

int main() {
    ComplexPlot data1;  // step
    ComplexPlot data2;  // step + noise
    ComplexPlot data3;  // sin1
    ComplexPlot data4;  // sin10 + sin1
    ComplexPlot data5;  // sin1 + noise
    ComplexPlot data6;  // lalf sin1, half sin3
    ComplexPlot data7;  // shift phase

    ComplexPlot transform1;
    ComplexPlot transform2;
    ComplexPlot transform3;
    ComplexPlot transform4;
    ComplexPlot transform5;
    ComplexPlot transform6;
    ComplexPlot transform7;
    ComplexPlot transformCheck;


    int quantity = 500;
    for(int i = 0; i < quantity; i++) {
        double x = (double)i/quantity * 15;

        if(i < quantity / 2 || i > 6 * quantity / 10)
            data1.push(x, {0, 0});
        else
            data1.push(x, {1, 0});

        data3.push(x, {sin(10*x), 0});
        data4.push(x, {sin(20*x) + sin(10*x), 0});

        if(i < quantity / 2)
            data6.push(x, {sin(9*x), 0});
        else
            data6.push(x, {sin(81*x), 0});

        if(i < quantity / 2)
            data7.push(x, {sin(5*x), 0});
        else
            data7.push(x, {sin(-5*x), 0});
    }
    data2 = data1;
    addNoise(data2);


    //addNoise(data, 1);
    fourierTransform(data1, transform1);
    fourierTransform(data2, transform2);
    fourierTransform(data3, transform3);
    fourierTransform(data4, transform4);
    fourierTransform(data6, transform6);
    fourierTransform(data7, transform7);
    dft1(data1, transformCheck);


    saveVectorPoint2DToFile(data1.real(), "data1.dat");
    saveVectorPoint2DToFile(transform1.abs(), "transform1.dat", transform1.size()/2);
    saveVectorPoint2DToFile(data2.real(), "data2.dat");
    saveVectorPoint2DToFile(transform2.abs(), "transform2.dat", transform2.size()/2);
    saveVectorPoint2DToFile(data3.real(), "data3.dat");
    saveVectorPoint2DToFile(transform3.abs(), "transform3.dat", transform3.size()/2);
    saveVectorPoint2DToFile(data4.real(), "data4.dat");
    saveVectorPoint2DToFile(transform4.abs(), "transform4.dat", transform4.size()/2);
    saveVectorPoint2DToFile(data6.real(), "data6.dat");
    saveVectorPoint2DToFile(transform6.abs(), "transform6.dat", transform6.size()/2);
    saveVectorPoint2DToFile(data7.real(), "data7.dat");
    saveVectorPoint2DToFile(transform7.abs(), "transform7.dat", transform7.size()/8);

    saveVectorPoint2DToFile(transform1.phase(), "arg1.dat");
    saveVectorPoint2DToFile(transform6.phase(), "arg2.dat");

    GnuplotPipe gp;
    //gp.sendLine(R"(set xrange [-10:10])");
    //gp.sendLine(R"(set yrange [-2:2])");
    gp.sendLine(R"(set multiplot layout 3,4)");
    gp.sendLine(R"(unset key)");

    gp.sendLine(R"(plot "data1.dat" with lines)");
    gp.sendLine(R"(plot "transform1.dat" with lines)");
    gp.sendLine(R"(plot "data2.dat" with lines)");
    gp.sendLine(R"(plot "transform2.dat" with lines)");
    //gp.sendLine(R"(plot "arg1.dat" with lines)");

    gp.sendLine(R"(plot "data3.dat" with lines)");
    gp.sendLine(R"(plot "transform3.dat" with lines)");
    gp.sendLine(R"(plot "data4.dat" with lines)");
    gp.sendLine(R"(plot "transform4.dat" with lines)");
    //gp.sendLine(R"(plot "arg2.dat" with lines)");

    gp.sendLine(R"(plot "data6.dat" with lines)");
    gp.sendLine(R"(plot "transform6.dat" with lines)");
    gp.sendLine(R"(plot "data7.dat" with lines)");
    gp.sendLine(R"(plot "transform7.dat" with lines)");

    return 0;
}