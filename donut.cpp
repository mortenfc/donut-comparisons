#include <utility>
#include <math.h>
#include <stdio.h>
#include <array>
#include <unistd.h>
#include <cmath>

constexpr uint32_t screen_width = 75;
constexpr uint32_t screen_height = 75;

constexpr float theta_spacing = 0.07;
constexpr float phi_spacing = 0.02;

constexpr float R1 = 3;
constexpr float R2 = 5;
constexpr float K2 = 20;
constexpr float size_scaler = 0.6;
// Calculate K1 based on screen size: the maximum x-distance occurs
// roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
// want that to be displaced 3/8ths of the width of the screen, which
// is 3/4th of the way from the center to the side of the screen.
// screen_width*3/8 = K1*(R1+R2)/(K2+0)
// screen_width*K2*3/(8*(R1+R2)) = K1
constexpr float K1 = screen_width * K2 * 3 / (8 * (R1 + R2)) * size_scaler;

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

static auto output{make_2d_array<char>(' ')};
static auto output_init_copy{output};

static auto zbuffer{make_2d_array<float>(0.0)};
static auto zbuffer_init_copy{zbuffer};

// Const because of viewer direction. Z buffer removes coordinates behind what's closest to the viewer.
constexpr float Lz = -1;

void render_frame(float A, float B, float Lx, float Ly)
{

    output = output_init_copy;
    zbuffer = zbuffer_init_copy;
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

            float L = Lz * (sinA * sintheta + cosA * costheta * sinphi) + Ly * (costheta * cosphi * sinB + cosB * (cosA * sintheta - costheta * sinA * sinphi)) + Lx * (cosB * costheta * cosphi - sinB * (cosA * sintheta - costheta * sinA * sinphi));
            // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
            // is pointing away from us, so we won't bother trying to plot it.
            if (L > -1)
            {
                L++;
                // test against the z-buffer.  larger 1/z means the pixel is
                // closer to the viewer than what's already plotted.
                if (ooz > zbuffer[xp][yp])
                {
                    zbuffer[xp][yp] = ooz;
                    int luminance_index = L * 2.2; // max index is 5
                    // luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
                    // now we lookup the character corresponding to the
                    // luminance and plot it in our output:
                    output[xp][yp] = ".-=?%@"[luminance_index];

                    output[light_source_x][light_source_y] = 'X';
                }
            }
        }
    }

    // now, dump output[] to the screen.
    // bring cursor to "home" location, in just about any currently-used
    // terminal emulation mode
    printf("\x1b[H");
    for (uint32_t j = 0; j < screen_height; j++)
    {
        for (uint32_t i = 0; i < screen_width; i++)
        {
            putchar(output[i][j]);
        }
        putchar('\n');
    }
}

int main(int argc, char const *argv[])
{
    float duration_s = 5;
    float steps = 3000;
    float sleep_dur = duration_s / steps;
    float inc = M_PI / steps;

    float i{0}, j{0}, Ly{0}, Lx{0}, i_ls{0};
    printf("\033c");
    while (true)
    {
        i_ls += inc * 0.95;
        i += inc * 1.55;
        j += inc * 0.55;
        Lx = cos(i_ls);
        Ly = sin(i_ls);
        render_frame(i, j, Lx, Ly);
        sleep(sleep_dur);
    }

    return 0;
}
