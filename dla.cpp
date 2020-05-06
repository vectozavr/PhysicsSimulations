//
// Created by ivan- on 03.05.2020.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include "vemath.h"

using namespace std;
using namespace vemath;

[[nodiscard]] Point2D randomVelocity2D(double abs) {
    Point2D randomDir;

    double angle = 2*PI*(double)rand()/RAND_MAX;

    randomDir.x = abs*cos(angle);
    randomDir.y = abs*sin(angle);

    return randomDir;
}

int main() {
    srand(126);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 4;

    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "turtles", sf::Style::Default, settings);
    window.clear(sf::Color(255, 255, 255, 0));
    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    unsigned long long N = 1;
    double R = 500*100;
    double velocity = 500;
    double acceleration = 1.f;

    int i_totalTime = 0;
    double totalTime = 0;

    int steps = 0;

    double scale = 0.01;

    vector<Point2D> turtles{N, {0,0}};

    turtles[0].x = 0;
    turtles[0].y =  0;

    vector<Point2D> history;

    while (window.isOpen())
    {
        vector<Point2D> lastPositions = turtles;

        tp2 = chrono::system_clock::now();
        chrono::duration <double> elapsedTime = tp2 - tp1;
        tp1 = tp2;
        double d_elapsedTime = elapsedTime.count();

        std::string title = "turtles " + std::to_string((double)1/elapsedTime.count());
        window.setTitle(title);
        sf::Event event{};
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        //cout << totalTime << endl;
        //window.clear();     // отчистка

        window.clear(sf::Color(255, 255, 255, 0));
        sf::RectangleShape rectangle(sf::Vector2f(2*R*SCALE, 2*R*SCALE));
        rectangle.move(SCREEN_WIDTH/2 - R*SCALE, SCREEN_HEIGHT/2 - R*SCALE);
        rectangle.setFillColor(sf::Color(255, 255, 255));
        rectangle.setOutlineThickness(10.f); // устанавливаем толщину контура круга
        rectangle.setOutlineColor(sf::Color(80,220,50)); // устанавливаем цвет контура
        //window.draw(rectangle);

        for(int i = 0; i < turtles.size(); i++) {

            int circleRadius = 1;

            sf::CircleShape circle(circleRadius);

            circle.setFillColor(sf::Color(255, 255, 255, 50));

            circle.setPosition(SCREEN_WIDTH/2 + (int)turtles[i].x*SCALE - circleRadius/2, SCREEN_HEIGHT/2 + (int)turtles[i].y*SCALE - circleRadius/2);
            //window.draw(circle);

            for(int s = steps; s < history.size()-1; s++) {
                sf::Vertex line[] =
                        {
                                sf::Vertex(sf::Vector2f(SCREEN_WIDTH / 2 + (int) history[s].x * scale,
                                                        SCREEN_HEIGHT / 2 + (int) history[s].y * scale)),
                                sf::Vertex(sf::Vector2f(SCREEN_WIDTH / 2 + (int) history[s+1].x * scale,
                                                        SCREEN_HEIGHT / 2 + (int) history[s+1].y * scale))
                        };

                line->color = sf::Color{0, 0, 0, 50};

                window.draw(line, 2, sf::Lines);
            }
            //steps = history.size()-1;

        }
        i_totalTime = (int) totalTime;
        window.display();   // отображение

    }

    system("pause");

    return 0;
}
