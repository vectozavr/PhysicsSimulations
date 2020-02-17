//
// Created by ivan- on 17.02.2020.
//

#include "mymath.h"

using namespace mymath;

bool mymath::saveVectorPoint2DToFile(const std::vector<mymath::Point2D> &data, const std::string &fileName) {
    std::ofstream _ofstream(fileName);
    if (!_ofstream.is_open())
        return false;

    for (auto d : data)
        _ofstream << d.x << "   " << d.y << std::endl;

    _ofstream.close();
    return true;
}

void mymath::fourierTransform(const ComplexPlot& data, ComplexPlot& transform) {
    transform.v_c.clear();
    for (int k = 0; k < data.v_c.size(); k++) {
        mymath::Complex transformed = {0, 0};
        for (int n = 0; n < data.v_c.size(); n++) {
            double p = 2 * PI * k * n / data.v_c.size();
            Complex tr = {cos(p) / data.v_c.size(), -sin(p) / data.v_c.size()};
            transformed += data.v_c[n].second * tr;
        }
        transform.push(2 * PI * k / data.v_c.size(), transformed);
    }
}

void mymath::inverseFourierTransform(const ComplexPlot& data, ComplexPlot& transform) {
    transform.v_c.clear();
    for (int k = 0; k < data.v_c.size(); k++) {
        mymath::Complex transformed = {0, 0};
        for (int n = 0; n < data.v_c.size(); n++) {
            double p = 2 * PI * k * n / data.v_c.size();
            Complex tr = {cos(p), sin(p)};
            transformed += data.v_c[n].second * tr;
        }
        transform.push(2 * PI * k / data.v_c.size(), transformed);
    }
}

void mymath::addNoise(ComplexPlot& data, double noiseAmplitude) {
    srand(1235);
    for(auto& d : data.v_c) {
        d.second.real += noiseAmplitude * rand()/RAND_MAX - noiseAmplitude/2;
        d.second.imagine += noiseAmplitude * rand()/RAND_MAX - noiseAmplitude/2;
    }
}
