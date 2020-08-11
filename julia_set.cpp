//
// Created by ivan- on 06.05.2020.
//

// LMB - move camera
// Mouse Wheel - scale
// Spacebar - increase iteration
// LeftShift+Spacebar - decrease iteration

#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <cmath>
#include <iostream>
#include <string>

using namespace sf;
using namespace std;

// General constatns
const int WIDTH = 800, HEIGHT = 800, MAX_IT = 200;
const double MAX_LEN = 2.0;

// General variables
static Vector2<double> offset = Vector2<double>(0, 0);
static Vector2i mousePos = Vector2i(0, 0);
static double scale = 0.25f;

// CPU rendering variables
static Vertex* allMap = nullptr;

// GPU rendering variables
static Shader* shader = nullptr;
static Texture* texture = nullptr;
static Sprite* screen = nullptr;

// Shader for OpenGL 4.0 and higher
static string double_shader_text =
"#version 400\n\
precision highp float;\n\
uniform float limit;\n\
uniform vec2 offset1;\n\
uniform vec2 offset2;\n\
uniform vec2 vscale;\n\
uniform vec2 res;\n\
uniform vec2 c;\n\
\n\
void main()\n\
{\n\
	const int maxit = 255;\n\
    dvec2 z = gl_FragCoord.xy / res;\n\
	dvec2 offset = offset1;\n\
	offset += offset2;\n\
	z = (z - dvec2(0.5, 0.5)) / vscale + offset;\n\
    int k = maxit;\n\
    for (int i = 0; i < maxit; i++)\n\
    {\n\
        z = dvec2(z.x * z.x - z.y * z.y + c.x, z.y * z.x * 2.0 + c.y);\n\
        if (z.x * z.x + z.y * z.y > limit)\n\
        {\n\
            k = i;\n\
            break;\n\
        }\n\
    }\n\
    float col = 1.0 - float(k) / float(maxit);\n\
    gl_FragColor = vec4(col, col, col, 1.0);\n\
}";

// Shader for OpenGL 2.0 and higher (but can only use float, so high accuracy losses)
static string float_shader_text =
"#version 110\n\
precision highp float;\n\
uniform float limit;\n\
uniform vec2 offset;\n\
uniform vec2 vscale;\n\
uniform vec2 res;\n\
uniform vec2 c;\n\
\n\
void main()\n\
{\n\
	const int maxit = 255;\n\
    vec2 z = gl_FragCoord.xy / res;\n\
	z = (z - vec2(0.5, 0.5)) / vscale + offset;\n\
    int k = maxit;\n\
    for (int i = 0; i < maxit; i++)\n\
    {\n\
        z = vec2(z.x * z.x - z.y * z.y + c.x, z.y * z.x * 2.0 + c.y);\n\
        if (z.x * z.x + z.y * z.y > limit)\n\
        {\n\
            k = i;\n\
            break;\n\
        }\n\
    }\n\
    float col = 1.0 - float(k) / float(maxit);\n\
    gl_FragColor = vec4(col, col, col, 1.0);\n\
}";

struct C_complex
{
    double re, im;

    C_complex(double real = 0, double imaginary = 0)
    {
        im = imaginary;
        re = real;
    }

    C_complex operator * (const C_complex& cur)
    {
        double _re = this->re * cur.re - this->im * cur.im;
        double _im = this->im * cur.re + this->re * cur.im;
        return C_complex(_re, _im);
    }

    C_complex operator + (const C_complex& cur)
    {
        return C_complex(this->re + cur.re, this->im + cur.im);
    }

    double get_len()
    {
        return sqrt(re * re + im * im);
    }
};

inline double get_x(int x)
{
    return (static_cast<double>(x) / WIDTH - 0.5) / scale + offset.x;
}

inline double get_y(int y)
{
    return (static_cast<double>(y) / HEIGHT - 0.5) / scale - offset.y;
}

void draw_fractal(RenderWindow& window, const C_complex& c_0)
{
    int cur = 0;
    for (int i = 0; i < WIDTH; ++i)
        for (int j = 0; j < HEIGHT; ++j)
        {
            int cnt = MAX_IT;
            double x = get_x(i);
            double y = get_y(j);
            C_complex z(x, y);

            for (int it = 0; it < MAX_IT; ++it)
            {
                z = z * z + c_0;
                if (z.get_len() >= MAX_LEN)
                {
                    cnt = it;
                    break;
                }
                else
                {
                    cnt = cnt;
                }
            }

            int c = 255 - 255 * cnt / MAX_IT;

            allMap[cur++].color = Color(c, c, c);
        }
    window.draw(allMap, cur, PrimitiveType::Points);
}

void draw_fractal_opengl(RenderWindow& window, const C_complex& c_0)
{
    Glsl::Vec2 vscale = Glsl::Vec2(scale * window.getSize().y / window.getSize().x, scale);

    shader->setUniform("res", Glsl::Vec2(window.getSize().x, window.getSize().y));
    shader->setUniform("c", Glsl::Vec2(c_0.re, c_0.im));
    shader->setUniform("vscale", vscale);
    shader->setUniform("offset", Glsl::Vec2(offset.x, offset.y));
    shader->setUniform("limit", 10000.0f);

    window.draw(*screen, shader);
}

void draw_fractal_opengl_4(RenderWindow& window, const C_complex& c_0)
{
    Glsl::Vec2 vscale = Glsl::Vec2(scale * window.getSize().y / window.getSize().x, scale);
    Glsl::Vec2 offset1 = Glsl::Vec2(static_cast<int>(offset.x * 1e7) / 1e7, static_cast<int>(offset.y * 1e7) / 1e7);
    Glsl::Vec2 offset2 = Glsl::Vec2(offset.x - offset1.x, offset.y - offset1.y);

    shader->setUniform("res", Glsl::Vec2(window.getSize().x, window.getSize().y));
    shader->setUniform("c", Glsl::Vec2(c_0.re, c_0.im));
    shader->setUniform("vscale", vscale);
    shader->setUniform("offset1", offset1);
    shader->setUniform("offset2", offset2);
    shader->setUniform("limit", 10000.0f);

    window.draw(*screen, shader);
}

int main(int argc, char** argv)
{
    ContextSettings settings;
    settings.antialiasingLevel = 8;
    RenderWindow window(VideoMode(WIDTH, HEIGHT), "SFML works!", Style::Default, settings);
    window.setFramerateLimit(60);

    bool accuracy = argc <= 1 || strcmp(argv[1], "-s");

    Clock clock;
    int accPerSec = 30;
    int iter = 3891000;

    void (*draw_function)(RenderWindow&, const C_complex&);

    // OpenGL 4.0 shader. Fast and use double
    if (window.getSettings().majorVersion >= 4)
    {
        draw_function = draw_fractal_opengl_4;

        double_shader_text.replace(double_shader_text.find("maxit = 255") + 8, 3, to_string(MAX_IT));

        shader = new Shader();
        shader->loadFromMemory(double_shader_text, Shader::Type::Fragment);

        texture = new Texture();
        texture->create(1, 1);

        screen = new Sprite(*texture, IntRect(0, 0, WIDTH, HEIGHT));
    }
    // OpenGL 2.0 shader. Faster, but use float
    else if (!accuracy && window.getSettings().majorVersion >= 2)
    {
        draw_function = draw_fractal_opengl;

        float_shader_text.replace(float_shader_text.find("maxit = 255") + 8, 3, to_string(MAX_IT));

        shader = new Shader();
        shader->loadFromMemory(float_shader_text, Shader::Type::Fragment);

        texture = new Texture();
        texture->create(1, 1);

        screen = new Sprite(*texture, IntRect(0, 0, WIDTH, HEIGHT));
    }
    // CPU calucations. Slow, but use double and don't need specific OpenGL version
    else
    {
        draw_function = draw_fractal;

        glPointSize(1.0f);

        allMap = new Vertex[WIDTH * HEIGHT];
        for (int i = 0; i < WIDTH; i++)
            for (int j = 0; j < HEIGHT; j++)
                allMap[i * WIDTH + j] = Vertex(Vector2f(i, j), Color(255, 255, 255));
    } 

    while (window.isOpen())
    {
        // Events processing (close and mouse wheel move)
        Event event;
        while (window.pollEvent(event))
        {
            if (event.type == Event::Closed)
            {
                window.close();
            }
            if (event.type == Event::EventType::MouseWheelMoved && window.hasFocus())
            {
                scale *= pow(1.1, event.mouseWheel.delta);

                if (scale < 0.01f)
                    scale = 0.01f;
            }
        }

        // Input processing
        if (window.hasFocus())
        {
            // Iteration change
            if (Keyboard::isKeyPressed(Keyboard::Key::Space))
            {
                if (Keyboard::isKeyPressed(Keyboard::Key::LShift))
                    iter -= static_cast<int>(clock.getElapsedTime().asSeconds() * accPerSec);
                else
                    iter += static_cast<int>(clock.getElapsedTime().asSeconds() * accPerSec);
            }
            else
            {
                clock.restart();
            }
            // Offset change
            if (Mouse::isButtonPressed(Mouse::Button::Left))
            {
                Vector2i diff = Mouse::getPosition() - mousePos;
                offset.x -= diff.x / 1.0 / scale / window.getSize().x;
                offset.y += diff.y / 1.0 / scale / window.getSize().y;
            }
        }
        mousePos = Mouse::getPosition();

        // Draw
        window.clear();

        double a = (double)iter / 100000;
        double b = a / 100;

        draw_function(window, { b * sin(a), b * cos(a) });
        
        window.display();
    }

    // Delete points array
    if (allMap)
        delete[] allMap;

    cout << iter << endl;
    return 0;
}
