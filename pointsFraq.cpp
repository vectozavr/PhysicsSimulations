//
// Created by ivan- on 09.05.2020.
//

#include <SFML/Graphics.hpp>
#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

struct Point2D {
    double x = 0;
    double y = 0;
};

int randomPoint(unsigned count) {
    return (int)count*rand()/RAND_MAX;
}

int whatTriangle(double probability1) {
    if(rand() > probability1*RAND_MAX) {
        return 2;
    } else
        return 1;
}

void nextPoint(double& x, double& y) {
    float nextX, nextY;
    float r = (float)rand()/RAND_MAX;
    if (r < 0.01) {         // Рост стебля
        nextX =  0;
        nextY =  0.16 * y;
    } else if (r < 0.86) {  // Рост маленьких листьев
        nextX =  0.85 * x + 0.04 * y;
        nextY = -0.0 * x + 0.85 * y + 1.6;
    } else if (r < 0.93) {  // Рост больших листьев справа
        nextX =  0.20 * x - 0.26 * y;
        nextY =  0.23 * x + 0.22 * y + 1.6;
    } else {                // Рост больших листьев слева
        nextX = -0.15 * x + 0.28 * y;
        nextY =  0.26 * x + 0.24 * y + 0.44;
    }
    x = nextX;
    y = nextY;
}

int main() {
    int SCREEN_WIDTH    = 1920;
    int SCREEN_HEIGHT   = 1080;

    srand(126);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 4;

    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "point_frac", sf::Style::Default, settings);
    window.clear(sf::Color(255, 255, 255, 0));
    auto startTime = chrono::system_clock::now();
    auto tp1 = chrono::system_clock::now();
    auto tp2 = chrono::system_clock::now();

    double r0 = 3;
    unsigned long long N    = 3;
    double totalTime        = 0;
    double scale            = 80;

    double probability1 = 0.8;

    int r = 1;
    int R = 8;

    vector<Point2D> tr1;
    vector<Point2D> tr2;

    tr1.push_back({-0.95, 1.18 - 3});
    tr1.push_back({-1.74, 2.95 - 3});
    tr1.push_back({-4.98, 4.05 - 3});

    tr2.push_back({0,  0 - 3});
    tr2.push_back({6, 4 - 3});
    tr2.push_back({7.5, 12.7 - 3});

    vector<Point2D> points;
    //points.push_back({0, 3});
    //points.push_back({-3.5, -1.5});
    //points.push_back({3, -3});
    //points.push_back({2, 1});
    //points.push_back({3, -2});
    //points.push_back({1.5, -3});
    for(int i = 0; i < N; i++) {
        points.push_back({r0*cos(i*2*3.1415/N), r0*sin(i*2*3.1415/N)});
    }

    Point2D cur = {(points[0].x + points[1].x)/2, (points[0].y + points[1].y)/2};

    //for(int k = 0; k < points.size(); k++) {
    //    sf::CircleShape circle(R);
    //    circle.setFillColor(sf::Color(255, 0, 0, 155));
    //    circle.setPosition(SCREEN_WIDTH/2 + points[k].x*scale - R, SCREEN_HEIGHT/2 - points[k].y*scale - R);
    //    window.draw(circle);
    //}

    sf::ConvexShape triangle1;
    sf::ConvexShape triangle2;
    triangle1.setPointCount(3);
    triangle2.setPointCount(3);
    for(int k = 0; k < 3; k++) {
        triangle1.setPoint(k, sf::Vector2f(tr1[k].x*scale, -tr1[k].y*scale));
        triangle2.setPoint(k, sf::Vector2f(tr2[k].x*scale, -tr2[k].y*scale));
    }
    triangle1.setFillColor({255, 0, 0, 100});
    triangle1.setPosition(SCREEN_WIDTH/2, SCREEN_HEIGHT/2);
    triangle2.setFillColor({0, 255, 0, 100});
    triangle2.setPosition(SCREEN_WIDTH/2, SCREEN_HEIGHT/2);
    //window.draw(triangle1);
    //window.draw(triangle2);


    int flag = 10000;

    while (window.isOpen())
    {
        tp2 = chrono::system_clock::now();
        chrono::duration <double> elapsedTime = tp2 - tp1;
        tp1 = tp2;
        double d_elapsedTime = elapsedTime.count();

        std::string title = "point_frac " + std::to_string((double)1/elapsedTime.count());
        window.setTitle(title);
        sf::Event event{};
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        //if(sf::Keyboard::isKeyPressed(sf::Keyboard::Space) && (flag > 1500))
        {
            flag = 0;
            sf::CircleShape circle(r);
            circle.setFillColor(sf::Color(255, 0, 0, 100));
            circle.setPosition(SCREEN_WIDTH / 2 + cur.x * scale,
                               SCREEN_HEIGHT / 2 - cur.y * scale + 350);
            window.draw(circle);

            int rand = randomPoint(tr1.size());
            Point2D vecTranslation;
            if(whatTriangle(probability1) == 1)
                vecTranslation = {(tr1[rand].x - cur.x) / 2, (tr1[rand].y - cur.y) / 2};
            else
                vecTranslation = {(tr2[rand].x - cur.x) / 2, (tr2[rand].y - cur.y) / 2};

            //cur.x += vecTranslation.x;
            //cur.y += vecTranslation.y;
            nextPoint(cur.x, cur.y);
        }
        window.display();   // отображение
        flag++;
        totalTime++;

    }

    system("pause");
    return 0;
}
