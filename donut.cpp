#include <utility>
#include <math.h>
#include <stdio.h>
#include <array>
#include <vector>
#include <unistd.h>
#include <unistd.h>
#include <iostream>
#include <stdint.h>

#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// TODO Debug why big values of these causes black dots to appear.
//* Theory: Does not calculate the same xyz for bigger resolutions / larger angle spacings,
//* so points between points will be calculated as rotating, meaning ooz can print a point behind a point.
constexpr uint32_t screen_width = 50;
constexpr uint32_t screen_height = 50;
constexpr uint8_t opengl_resolution_mutliplier = 7;

constexpr float theta_spacing = 0.07;
constexpr float phi_spacing = 0.02;

constexpr float R1 = 3;
constexpr float R2 = 5;
constexpr float K2 = 20;
constexpr float size_scaler = 1;
// Calculate K1 based on screen size: the maximum x-distance occurs
// roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
// want that to be displaced 3/8ths of the width of the screen, which
// is 3/4th of the way from the center to the side of the screen.
// screen_width*3/8 = K1*(R1+R2)/(K2+0)
// screen_width*K2*3/(8*(R1+R2)) = K1
constexpr float K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) * size_scaler;

// Const because of viewer direction. Z buffer removes coordinates behind what's closest to the viewer.
constexpr float Lz = -0.75;

static GLFWwindow *g_window;

static void error_callback(int error, const char *description)
{
    std::cerr << "Error: " << description << std::endl;
}

bool gl_prepare()
{
    glfwInit();
    glfwSetErrorCallback(error_callback);
    /* Create a windowed mode g_window and its OpenGL context */
    g_window = glfwCreateWindow(screen_width * opengl_resolution_mutliplier, screen_height * opengl_resolution_mutliplier, "Donut OpenGL", NULL, NULL);
    if (!g_window)
    {
        std::cout << "No window created ... " << std::endl;
        glfwTerminate();
        return false;
    }

    /* Make the g_window's context current */
    glfwMakeContextCurrent(g_window);

    return true;
}

template <typename T>
constexpr auto make_2d_array(T init)
{
    std::array<std::array<T, screen_height>, screen_width> matrix;

    for (int i = 0; i < screen_width; i++)
    {
        for (int j = 0; j < screen_height; j++)
        {
            matrix[i][j] = init;
        }
    }

    return matrix;
}

static auto ascii{make_2d_array<char>(' ')};
static auto ascii_init_copy{ascii};

static auto zbuffer{make_2d_array<float>(0.0)};
static auto zbuffer_init_copy{zbuffer};

static auto lum{make_2d_array<float>(0.0)};
static auto lum_init_copy{lum};

void drawPolygon(float x, float y, float radius, int num_corners)
{
    // Requires GL_TRIANGLES
    glBegin(GL_TRIANGLES);
    for (int corner = 0; corner < num_corners; corner++)
    {
        float theta1 = 2.0f * M_PI * float(corner) / float(num_corners);
        float theta2 = 2.0f * M_PI * float(corner + 1) / float(num_corners);
        glVertex2f(x, y);
        glVertex2f(x + radius * cosf(theta1), y + radius * sinf(theta1));
        glVertex2f(x + radius * cosf(theta2), y + radius * sinf(theta2));
    }
    glEnd();
}

using Ascii = std::array<std::array<char, screen_height>, screen_width>;
using Lums = std::array<std::array<float, screen_height>, screen_width>;
void glDisplayPoints(
    Lums const &points)
{
    // https://eng.libretexts.org/Bookshelves/Computer_Science/Applied_Programming/Book%3A_Introduction_to_Computer_Graphics_(Eck)/03%3A_OpenGL_1.1-_Geometry/3.01%3A_Shapes_and_Colors_in_OpenGL_1.1

    constexpr static float sqrt_2 = sqrtf(2.0F);

    glClear(GL_COLOR_BUFFER_BIT);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(opengl_resolution_mutliplier);

    glBegin(GL_POINTS);
    for (float j = 0; j < float(screen_height); j++)
    {
        for (float i = 0; i < float(screen_width); i++)
        {
            float x = 2.0f * i / screen_width - 1.0f;
            float y = 2.0f * j / screen_height - 1.0f;

            if (points[i][j] == -10)
            {
                glColor3f(1, 1, 0); // Yellow
                constexpr static int num_corners = 20;
                constexpr static float radius = 0.002 * (5 + opengl_resolution_mutliplier);
                glEnd();
                drawPolygon(x, y, radius, num_corners);
                glBegin(GL_POINTS);
                continue;
            }
            const float lum = points[i][j] / sqrt_2;
            if (lum > 0.F)
            {
                glColor3f(lum, lum, lum); // Gray
                glVertex2f(x, y);
            }
        }
    }
    glEnd();

    /* Swap front and back buffers */
    glfwSwapBuffers(g_window);

    /* Poll for and process events */
    glfwPollEvents();
}

void render_frame(float A, float B, float Lx, float Ly)
{
    ascii = ascii_init_copy;
    zbuffer = zbuffer_init_copy;
    lum = lum_init_copy;
    // precompute sines and cosines of A and B
    float cosA = cos(A), sinA = sin(A);
    float cosB = cos(B), sinB = sin(B);

    // theta goes around the cross-sectional circle of a torus
    for (float theta = 0; theta < 2.0F * M_PI; theta += theta_spacing)
    {
        // precompute sines and cosines of theta
        float costheta = cos(theta), sintheta = sin(theta);

        // phi goes around the center of revolution of a torus
        for (float phi = 0; phi < 2 * M_PI; phi += phi_spacing)
        {
            // precompute sines and cosines of phi
            float cosphi = cos(phi), sinphi = sin(phi);

            // the x,y coordinate of the circle, before revolving (factored
            // out of the above equations)
            float circlex = R2 + R1 * costheta;
            float circley = R1 * sintheta;

            // final 3D (x,y,z) coordinate after rotations, directly from
            // our math above
            float x = circlex * (cosB * cosphi + sinA * sinB * sinphi) - circley * cosA * sinB;
            float y = circlex * (sinB * cosphi - sinA * cosB * sinphi) + circley * cosA * cosB;
            float z = K2 + cosA * circlex * sinphi + circley * sinA;
            float ooz = 1.0F / z; // "one over z"

            // x and y projection.  note that y is negated here, because y
            // goes up in 3D space but down on 2D displays.
            int xp = (int)(screen_width / 2.0F + K1 * ooz * x);
            int yp = (int)(screen_height / 2.0F - K1 * ooz * y);

            int light_source_x = (int)(((float)screen_width / 2.0F) + Lx * ((float)screen_width / 2.01F));
            int light_source_y = (int)(((float)screen_height / 2.0F) - Ly * ((float)screen_height / 2.01F));

            const float L = Lz * (sinA * sintheta + cosA * costheta * sinphi) + Ly * (costheta * cosphi * sinB + cosB * (cosA * sintheta - costheta * sinA * sinphi)) + Lx * (cosB * costheta * cosphi - sinB * (cosA * sintheta - costheta * sinA * sinphi));
            // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
            // is pointing away from us, so we won't bother trying to plot it.
            if (L > -1)
            {
                // test against the z-buffer. larger 1/z means the pixel is
                // closer to the viewer than what's already plotted.
                if (ooz > zbuffer[xp][yp])
                {
                    zbuffer[xp][yp] = ooz;
                    int luminance_index = (L + 1.F) * 2.2; // max index is 5
                    // luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
                    // now we lookup the character corresponding to the
                    // luminance and plot it in our ascii:
                    ascii[xp][yp] = ".-=?%@"[luminance_index];
                    ascii[light_source_x][light_source_y] = 'X';

                    lum[xp][yp] = L;
                    lum[light_source_x][light_source_y] = -10;
                }
            }
        }
    }
}

void asciiDisplayPoints(Ascii const &ascii)
{
    // now, dump ascii[] to the screen.
    // bring cursor to "home" location, in just about any currently-used
    // terminal emulation mode
    printf("\x1b[H");
    for (uint32_t j = 0; j < screen_height; j++)
    {
        for (uint32_t i = 0; i < screen_width; i++)
        {
            putchar(ascii[i][j]);
        }
        putchar('\n');
    }
}

int main(int argc, char const *argv[])
{
    // return gl_main(argc, (char**)argv);
    if (!gl_prepare())
    {
        return -1;
    }
    float duration_s = 5;
    float steps = 3000;
    float sleep_dur = duration_s / steps;
    float inc = M_PI / steps;

    float i{0}, j{0}, Ly{0}, Lx{0}, i_ls{0};
    printf("\033c");
    while (!glfwWindowShouldClose(g_window))
    {
        i_ls += inc * 0.95;
        i += inc * 1.55;
        j += inc * 0.55;
        Lx = cos(i_ls);
        Ly = sin(i_ls);
        render_frame(i, j, Lx, Ly);
        // asciiDisplayPoints(ascii);
        glDisplayPoints(lum);
        usleep(int(sleep_dur) * 1000000);
    }

    glfwTerminate();
    return 0;
}
