//
// Created by ivan-Vectozavr on 17.02.2020.
//

#include "vemath.h"

using namespace vemath;

bool vemath::saveVectorPoint2DToFile(const std::vector<vemath::Point2D> &data, const std::string &fileName, unsigned long long N) {
    std::ofstream _ofstream(fileName);
    if (!_ofstream.is_open())
        return false;

    N = N==0 ? data.size() : N;
    for(int i = 0; i < N; i++) {
        _ofstream << data[i].x << "   " << data[i].y << std::endl;
    }

    _ofstream.close();
    return true;
}

void vemath::fourierTransform(const ComplexPlot& data, ComplexPlot& transform) {
    transform.v_c.clear();
    for (int k = 0; k < data.v_c.size(); k++) {
        std::complex<double> transformed = {0, 0};
        for (int n = 0; n < data.v_c.size(); n++) {
            double p = 2 * PI * k * n / data.v_c.size();
            std::complex<double> tr = {cos(p) / data.v_c.size(), -sin(p) / data.v_c.size()};
            transformed += data.v_c[n].second * tr;
        }
        transform.push(k, transformed);
    }
    transform._xn = data.v_c.back().first;
}

void vemath::inverseFourierTransform(const ComplexPlot& data, ComplexPlot& transform) {
    transform.v_c.clear();
    for (int k = 0; k < data.v_c.size(); k++) {
        std::complex<double> transformed = {0, 0};
        for (int n = 0; n < data.v_c.size(); n++) {
            double p = 2 * PI * k * n / data.v_c.size();
            std::complex<double> tr = {cos(p), sin(p)};
            transformed += data.v_c[n].second * tr;
        }
        transform.push(k*data.xn()/data.size(), transformed);
    }
}

void vemath::addNoise(ComplexPlot& data, double noiseAmplitude, int seed) {
    srand(seed);
    for(auto& d : data.v_c) {
        d.second += noiseAmplitude * rand()/RAND_MAX - noiseAmplitude/2;
    }
}

[[nodiscard]] Point3D randomDirection(int seed) {
    srand(seed);
    Point3D dir = {1, 1, 1};
    while(dir.abs() > 1) {
        dir.x = -1 + 2*rand()/RAND_MAX;
        dir.y = -1 + 2*rand()/RAND_MAX;
        dir.z = -1 + 2*rand()/RAND_MAX;
    }
    dir = dir.normalize();
    return dir;
}

//return min x from ComplexPlot
double vemath::minx(const ComplexPlot& data) {
    double min = data.v_c[0].first;
    for(auto comp : data.v_c)
        if(min > comp.first)
            min = comp.first;
    return min;
}
//return max x from ComplexPlot
double vemath::maxx(const ComplexPlot& data) {
    double max = data.v_c[0].first;
    for(auto comp : data.v_c)
        if(max < comp.first)
            max = comp.first;
    return max;
}
//return min y from ComplexPlot
double vemath::miny(const ComplexPlot& data) {
    double min = data.v_c[0].second.real();
    for(auto comp : data.v_c)
        if(min > comp.second.real())
            min = comp.second.real();
    return min;
}
//return max y from ComplexPlot
double vemath::maxy(const ComplexPlot& data) {
    double max = data.v_c[0].second.real();
    for(auto comp : data.v_c)
        if(max < comp.second.real())
            max = comp.second.real();
    return max;
}

void vemath::convolution(const ComplexPlot& data1, const ComplexPlot& data2, ComplexPlot& conv) {
    double min = std::min(minx(data1), minx(data2));
    double max = std::max(maxx(data1), maxx(data2));

    conv.v_c.clear();
    int N = std::max(data1.size(), data2.size());
    int n = std::min(data1.size(), data2.size());
    for(int i = 0; i < n + N; i++) {
        std::complex<double> result = {0, 0};
        for(int k = 0; k < data1.size() + data2.size(); k++) {
            if((k > data1.size()) || (i - k < 0) || (i - k > data2.size()))
                continue;
            result += data1.v_c[k].second * data2.v_c[i - k].second;
        }
        conv.push(i, result);
    }
}

void vemath::crossCorrelation(const ComplexPlot& data1, const ComplexPlot& data2, ComplexPlot& cross) {
    ComplexPlot inv_data1;
    for(int k = data1.size()-1; k >= 0; k--)
        inv_data1.push(data1.v_c[k].first, data1.v_c[k].second);

    convolution(inv_data1, data2, cross);
}

void vemath::highFilter(const ComplexPlot& in, ComplexPlot& out, int freq) {

}