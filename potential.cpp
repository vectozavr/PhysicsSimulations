//
// Created by ivan- on 24.02.2020.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include "vemath.h"
#include "settings.h"

using namespace std;
using namespace vemath;

void drawScene(sf::RenderWindow& window, vector<vector<sf::Color>>& scene) {
    for(int i = 0; i < scene.size(); i++) {
        for(int k = 0; k < scene[i].size(); k++) {
            sf::RectangleShape rect(sf::Vector2f(SCREEN_WIDTH / scene.size(), SCREEN_WIDTH / scene.size()));

            rect.setFillColor(scene[i][k]);
            rect.setPosition(i*SCREEN_WIDTH / scene.size(), k*SCREEN_WIDTH / scene.size());
            window.draw(rect);
        }
    }
}

double potentialCharges(Point3D position, vector<pair<Point3D, double>>& charges) {
    // Charges
    double potential = 0;
    for (auto ch : charges) {
        Point3D diff = position - ch.first;
        potential += ch.second / diff.abs();
    }
    return potential;
}

double potentialLines(Point3D position) {
    // Lines
    Point3D line1 = {position.x, 50, 0};
    Point3D line2 = {position.x, -50, 0};

    Point3D diff1 = position - line1;
    Point3D diff2 = position - line2;

    double potential = -log(diff1.abs()/150) - log(diff2.abs()/150);

    return potential;
}

double potentialTape(Point3D position) {
    // Tape
    double potential = 0;

    for(int i = -50; i < 50; i++) {
        Point3D line_i = {position.x, (double)i/100, 0};
        Point3D diff_i = position - line_i;
        potential += log(diff_i.abs()/100)/100;
    }

    return potential;
}

double potentialTapeAndLine(Point3D position) {
    // Tape
    double potential = 0;

    for(int i = -50; i < 50; i++) {
        Point3D line_i = {position.x, (double)i/20, 0};
        Point3D diff_i = position - line_i;
        potential += log(diff_i.abs()/100)/500;
    }
    Point3D line_1 = {position.x, -150, 0};
    Point3D diff_1 = position - line_1;
    potential -= log(diff_1.abs()/100)/10;

    return potential;
}

double twoPotentialTape(Point3D position) {
    // Tape
    double potential = 0;

    for(int i = -50; i < 50; i++) {
        Point3D line_i1 = {position.x, 50 + (double)i/20, 0};
        Point3D diff_i1 = position - line_i1;
        potential += log(diff_i1.abs()/100)/500;

        Point3D line_i2 = {position.x, -50 + (double)i/20, 0};
        Point3D diff_i2 = position - line_i2;
        potential += log(diff_i2.abs()/100)/500;
    }

    return potential;
}

int main() {
    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "potential");

    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    vector<vector<sf::Color>> scene;
    vector<pair<Point3D, double>> charges;

    int pixelSize = 3;
    double scale = (double)1;

    charges.push_back({{50, 0}, 1});
    charges.push_back({{-50, 0}, -1});

    bool redraw = true;

    while (window.isOpen()) {
        tp2 = chrono::system_clock::now();
        chrono::duration<double> elapsedTime = tp2 - tp1;
        tp1 = tp2;
        double d_elapsedTime = elapsedTime.count();

        std::string title = "potential " + std::to_string((double) 1 / elapsedTime.count());
        window.setTitle(title);
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        if(redraw) {
            for (int i = 0; i < SCREEN_WIDTH/pixelSize; i++) {
                vector<sf::Color> colors;
                for (int k = 0; k < SCREEN_HEIGHT/pixelSize; k++) {
                    Point3D position = {(i * pixelSize - (double) SCREEN_WIDTH / 2) * scale,(k * pixelSize - (double) SCREEN_HEIGHT / 2) * scale};

                    //double potential = potentialCharges(position, charges); // charges
                    //double potential = potentialLines(position);  // lines
                    //double potential = potentialTape(position);  // tape
                    //double potential = potentialTapeAndLine(position);  // tape and line
                    double potential = twoPotentialTape(position);  // two tapes

                    if(potential < 0) potential *= -1;
                    vector<int> res_color = colorInterpolate((potential*1));
                    colors.emplace_back(static_cast<sf::Uint8>(res_color[0]), static_cast<sf::Uint8>(res_color[1]), static_cast<sf::Uint8>(res_color[2]));
                }
                scene.push_back(colors);
            }
            redraw = false;
        }
        window.clear();     // отчистка
        drawScene(window, scene);
        window.display();   // отображение
    }

    return 0;
}
