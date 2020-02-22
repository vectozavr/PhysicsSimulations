//
// Created by ivan- on 16.02.2020.
//
#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include "vemath.h"

using namespace std;
using namespace vemath;

[[nodiscard]] Point3D randomDirection1(int seed) {
    Point3D dir = {1, 1, 1};
    while(dir.abs() > 1) {
        dir.x = -1 + 2*(double)rand()/RAND_MAX;
        dir.y = -1 + 2*(double)rand()/RAND_MAX;
        dir.z = -1 + 2*(double)rand()/RAND_MAX;
    }
    dir = dir.normalize();
    return dir;
}

[[nodiscard]] Point3D randomVelocity(double abs) {
    Point3D randomDir = randomDirection1(1);
    randomDir *= abs;
    return randomDir;
}

struct Rocket {
    double x = 0;
    double y = 0;
    double z = 0;

    double Vx = 0;
    double Vy = 0;
    double Vz = 0;

    double explosionHeight = 200;
    double explosionPower = 100;

    double mass = 0;
    int numberOfParticles = 400;

    double Ax = 0;
    double Ay = -9.81;
    double Az = 0;

    //sf::Color C_color = {dynamic_cast<sf::Utf8>(255*rand()/RAND_MAX), 255*rand()/RAND_MAX, 255*rand()/RAND_MAX};

    vector<Point3D> v_particles;

    void explode(double power = 100) {
        v_particles.clear();
        for(int i = 0; i < numberOfParticles; i++) {
            Point3D r_v = randomVelocity(power);
            v_particles.push_back({x, y, z, Vx + r_v.x, Vy + r_v.y, Vz + r_v.z, Ax, Ay, Az, mass/numberOfParticles});
        }
    }

    void draw(sf::RenderWindow& window, sf::Color color = {255, 255, 255}) {
        if(v_particles.empty()) {
            sf::CircleShape circle(2.f);
            circle.setFillColor(color);
            circle.setPosition(SCREEN_WIDTH/2 + x*SCALE - 1, SCREEN_HEIGHT - y*SCALE - 1);
            window.draw(circle);
        } else {
            for(auto p : v_particles) {
                sf::CircleShape circle(2.f);
                circle.setFillColor(color);
                circle.setPosition(SCREEN_WIDTH/2 + p.x*SCALE - 1, SCREEN_HEIGHT - p.y*SCALE - 1);
                window.draw(circle);
            }
        }
    }
};

double calculateShift(double elapsedTime, Rocket& rocket, double dt = 0.001) {
    double time = 0;
    while (time < elapsedTime) {
        if(rocket.v_particles.empty()) {
            //one whole rocket
            rocket.x += rocket.Vx*dt;
            rocket.y += rocket.Vy*dt;
            rocket.z += rocket.Vz*dt;
            rocket.Vx += rocket.Ax*dt;
            rocket.Vy += rocket.Ay*dt;
            rocket.Vz += rocket.Az*dt;
            if(rocket.y >= rocket.explosionHeight)
                rocket.explode(rocket.explosionPower);
        } else {
            for(auto& p : rocket.v_particles) {
                p.x     += p.Vx*dt;
                p.y     += p.Vy*dt;
                p.z     += p.Vz*dt;
                p.Vx    += p.Ax*dt;
                p.Vy    += p.Ay*dt;
                p.Vz    += p.Az*dt;
            }
        }
        time += dt;
    }
    return time;
}



int main() {
    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "salute");

    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    double rocketMass = 12*1000;
    double rocketVelocity = 100;
    double explosionHeight = 200;
    double explosionPower = 100;
    int numberOfParticles = 1000;

    int numberOfRockets = 5;

    vector<Rocket> v_rockets;
    for(int i = 0; i < numberOfRockets; i++)
        v_rockets.push_back({(double)SCREEN_WIDTH/SCALE*((double)(i+0.5)/numberOfRockets - 0.5), 0, 0, 0, rocketVelocity, 0, explosionHeight, explosionPower, rocketMass, numberOfParticles});

    //    double x = 0;
    //    double y = 0;
    //    double z = 0;
    //
    //    double Vx = 0;
    //    double Vy = 0;
    //    double Vz = 0;
    //
    //    double explosionHeight = 200;
    //    double explosionPower = 100;
    //
    //    double mass = 0;
    //    int numberOfParticles = 400;
    //
    //    double Ax = 0;
    //    double Ay = -9.81;
    //    double Az = 0;

    double totalTime = 0;
    while (window.isOpen())
    {
        tp2 = chrono::system_clock::now();
        chrono::duration <double> elapsedTime = tp2 - tp1;
        tp1 = tp2;
        double d_elapsedTime = elapsedTime.count();

        std::string title = "salute " + std::to_string((double)1/elapsedTime.count());
        window.setTitle(title);
        sf::Event event{};
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        double plus = 0;
        for(auto& r : v_rockets) {
            plus = calculateShift(d_elapsedTime, r);
            r.draw(window);
        }
        totalTime += plus;

        window.display();   // отображение
    }

    return 0;
}