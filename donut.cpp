#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <array>
#include <iostream>
#include <numeric>
#include <utility>

#include <chrono>
#include <functional>
#include <thread>

#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

//* Why do big resolutions + big angle spacing causes black dots to appear?
//* When projecting the 3d torus to 2d as XY, the rotational points won't overlap if it does not compute a point frequently enough,
//* leading to ooz can select a point on another surface between points in the surface facing the viewer.

// WARNING: Changing these rendering constants can heavily affect compilation duration

// One-frame donut rendering
constexpr uint32_t screen_width = 500;
constexpr uint32_t screen_height = screen_width;
constexpr uint8_t opengl_resolution_mutliplier = 1080 / screen_height;
constexpr float theta_spacing = 0.07 * (50.0 / screen_width);
constexpr float phi_spacing = 0.02 * (50.0 / screen_width);

// Rendering multiple frames (rotating the donut)
constexpr int kA_ls_den = 4;
constexpr int kA_den = 2;
constexpr int kB_den = 6;
constexpr int lcm(int a, int b, int c) {
  return std::lcm(std::lcm(a, b), c);
}
constexpr int rev_count = lcm(kA_ls_den, kA_den, kB_den);  // 12
constexpr float k(int denominator) {
  return 2.0f / denominator;
}
constexpr float kA_ls = k(kA_ls_den);  // 0.5;
constexpr float kA = k(kA_den);        // 1;
constexpr float kB = k(kB_den);        // 0.33333;
constexpr int its = 100 * (rev_count / 2);
constexpr float radians_one_cycle = float(rev_count) * M_PI;
constexpr float inc = radians_one_cycle / float(its);

// Torus dimensions
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

static GLFWwindow* g_window;

static void error_callback(int error, const char* description) {
  std::cerr << "Error: " << description << std::endl;
}

bool gl_prepare() {
  glfwInit();
  glfwSetErrorCallback(error_callback);
  /* Create a windowed mode g_window and its OpenGL context */
  g_window = glfwCreateWindow(screen_width * opengl_resolution_mutliplier, screen_height * opengl_resolution_mutliplier,
                              "Donut OpenGL", NULL, NULL);
  if (!g_window) {
    std::cout << "No window created ... " << std::endl;
    glfwTerminate();
    return false;
  }

  /* Make the g_window's context current */
  glfwMakeContextCurrent(g_window);

  return true;
}

template <typename T>
constexpr auto make_2d_array(T init) {
  std::array<std::array<T, screen_height>, screen_width> matrix;

  for (int i = 0; i < screen_width; i++) {
    for (int j = 0; j < screen_height; j++) {
      matrix[i][j] = init;
    }
  }

  return matrix;
}

void drawPolygon(float x, float y, float radius, int num_corners) {
  // Requires GL_TRIANGLES
  glBegin(GL_TRIANGLES);
  for (int corner = 0; corner < num_corners; corner++) {
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
void glDisplayPoints(Lums const& points) {
  // https://eng.libretexts.org/Bookshelves/Computer_Science/Applied_Programming/Book%3A_Introduction_to_Computer_Graphics_(Eck)/03%3A_OpenGL_1.1-_Geometry/3.01%3A_Shapes_and_Colors_in_OpenGL_1.1

  constexpr static float sqrt_2 = sqrtf(2.0F);

  glClear(GL_COLOR_BUFFER_BIT);

  glEnable(GL_POINT_SMOOTH);
  glPointSize(opengl_resolution_mutliplier);

  glBegin(GL_POINTS);
  for (float j = 0; j < float(screen_height); j++) {
    for (float i = 0; i < float(screen_width); i++) {
      float x = 2.0f * i / screen_width - 1.0f;
      float y = 2.0f * j / screen_height - 1.0f;

      if (points[i][j] == -10) {
        glColor3f(1, 1, 0);  // Yellow
        constexpr static int num_corners = 20;
        constexpr static float radius = 0.002 * (5 + opengl_resolution_mutliplier);
        glEnd();
        drawPolygon(x, y, radius, num_corners);
        glBegin(GL_POINTS);
        continue;
      }
      const float lum = points[i][j] / sqrt_2;
      if (lum > 0.F) {
        glColor3f(lum, lum, lum);  // Gray
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

constexpr auto ascii_init_copy{make_2d_array<char>(' ')};

constexpr auto zbuffer_init_copy{make_2d_array<float>(0.0)};

constexpr auto lum_init_copy{make_2d_array<float>(0.0)};

using DonutFrame = std::pair<Ascii, Lums>;

DonutFrame render_frame(float A, float B, float Lx, float Ly) {
  // TODO This math can be heavily optimized
  auto ascii = ascii_init_copy;
  auto zbuffer = zbuffer_init_copy;
  auto lum = lum_init_copy;
  // precompute sines and cosines of A and B
  float cosA = cos(A), sinA = sin(A);
  float cosB = cos(B), sinB = sin(B);

  // theta goes around the cross-sectional circle of a torus
  for (float theta = 0; theta < 2.0F * M_PI; theta += theta_spacing) {
    // precompute sines and cosines of theta
    float costheta = cos(theta), sintheta = sin(theta);

    // phi goes around the center of revolution of a torus
    for (float phi = 0; phi < 2 * M_PI; phi += phi_spacing) {
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
      float ooz = 1.0F / z;  // "one over z"

      // x and y projection.  note that y is negated here, because y
      // goes up in 3D space but down on 2D displays.
      int xp = (int)(screen_width / 2.0F + K1 * ooz * x);
      int yp = (int)(screen_height / 2.0F - K1 * ooz * y);

      int light_source_x = (int)(((float)screen_width / 2.0F) + Lx * ((float)screen_width / 2.01F));
      int light_source_y = (int)(((float)screen_height / 2.0F) - Ly * ((float)screen_height / 2.01F));

      const float L = Lz * (sinA * sintheta + cosA * costheta * sinphi) +
                      Ly * (costheta * cosphi * sinB + cosB * (cosA * sintheta - costheta * sinA * sinphi)) +
                      Lx * (cosB * costheta * cosphi - sinB * (cosA * sintheta - costheta * sinA * sinphi));
      // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
      // is pointing away from us, so we won't bother trying to plot it.
      if (L > -1) {
        // test against the z-buffer. larger 1/z means the pixel is
        // closer to the viewer than what's already plotted.
        if (ooz > zbuffer[xp][yp]) {
          zbuffer[xp][yp] = ooz;
          int luminance_index = (L + 1.F) * 2.2;  // max index is 5
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

  return std::make_pair(ascii, lum);
}

using DonutStorage = std::array<DonutFrame, its>;
constexpr DonutStorage render_donut_compile_time() {
  DonutStorage out;

  float A{0}, A_ls{0}, B{0};
  float i{0};
  for (float inc_cum{0.0}; inc_cum < radians_one_cycle; inc_cum += inc) {
    A_ls += inc * kA_ls;
    A += inc * kA;
    B += inc * kB;
    float Lx = cos(A_ls);
    float Ly = sin(A_ls);
    out[i] = render_frame(A, B, Lx, Ly);
    i++;
  }

  return out;
}

static DonutStorage asciis_and_lums{render_donut_compile_time()};

void asciiDisplayPoints(Ascii const& ascii) {
  // now, dump ascii[] to the screen.
  // bring cursor to "home" location, in just about any currently-used
  // terminal emulation mode
  printf("\x1b[H");
  for (uint32_t j = 0; j < screen_height; j++) {
    for (uint32_t i = 0; i < screen_width; i++) {
      putchar(ascii[i][j]);
    }
    putchar('\n');
  }
}

void sleep_for_frequency_Hz(float frequency, std::chrono::time_point<std::chrono::high_resolution_clock> start_time) {
  auto target_cycle_time = std::chrono::duration<double>(1.0 / frequency);
  auto processing_time = std::chrono::high_resolution_clock::now() - start_time;
  auto sleep_time_micros = std::chrono::duration_cast<std::chrono::microseconds>(target_cycle_time - processing_time);

  std::cout << "\nFrequency: " << frequency << "\nProcessing time in micros: "
            << std::chrono::duration_cast<std::chrono::microseconds>(processing_time).count()
            << "\nSleep time in micros: "
            << std::chrono::duration_cast<std::chrono::microseconds>(sleep_time_micros).count() << std::endl;
  // Sleep for the computed time
  if (sleep_time_micros.count() > 0) {
    std::this_thread::sleep_for(sleep_time_micros);
  }
  std::cout << "\033[2J\033[1;1H";
}

struct Args {
  enum class ProcessingMode { compile_time, run_time };
  ProcessingMode processing_mode = ProcessingMode::compile_time;

  enum class DisplayMode { ascii, opengl };
  DisplayMode display_mode = DisplayMode::opengl;
};

Args parse_arguments(int argc, char* const argv[]) {
  Args args;
  int opt;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      std::cout << "Options:\n"
                << "--display=ascii: Set display mode to ascii\n"
                << "--display=opengl: Set display mode to opengl (default)\n"
                << "--compute=runtime: Set compute mode to run time\n"
                << "--compute=compiletime: Set compute mode to compile time (default)\n";
      exit(0);
    } else if (strncmp(argv[i], "--display=", 10) == 0) {
      if (strcmp(argv[i] + 10, "ascii") == 0) {
        args.display_mode = Args::DisplayMode::ascii;
      } else if (strcmp(argv[i] + 10, "opengl") == 0) {
        args.display_mode = Args::DisplayMode::opengl;
      }
    } else if (strncmp(argv[i], "--compute=", 10) == 0) {
      if (strcmp(argv[i] + 10, "runtime") == 0) {
        args.processing_mode = Args::ProcessingMode::run_time;
      } else if (strcmp(argv[i] + 10, "compiletime") == 0) {
        args.processing_mode = Args::ProcessingMode::compile_time;
      }
    }
  }

  return args;
}

int main(int argc, char* const argv[]) {
  std::cout << "rev_count: " << rev_count << ", kA_ls: " << kA_ls << ", kA: " << kA << ", kB: " << kB
            << ", its: " << its << std::endl;
  Args args = parse_arguments(argc, argv);
  std::function<bool(int, DonutFrame const&)> display_func;
  std::function<DonutFrame const&(float, int)> compute_donut;

  switch (args.processing_mode) {
    case Args::ProcessingMode::run_time:
      compute_donut = [](float inc_sum, int i) -> DonutFrame const& {
        float A_ls = inc_sum * kA_ls;
        float A = inc_sum * kA;
        float B = inc_sum * kB;
        float Lx = cos(A_ls);
        float Ly = sin(A_ls);
        asciis_and_lums[i] = render_frame(A, B, Lx, Ly);
        return asciis_and_lums[i];
      };
      break;
    case Args::ProcessingMode::compile_time:
      compute_donut = [](float, int i) -> DonutFrame const& {
        return asciis_and_lums[i];
      };
      break;
  }

  switch (args.display_mode) {
    case Args::DisplayMode::ascii:
      printf("\033c");
      display_func = [](int i, DonutFrame const& donut) {
        asciiDisplayPoints(donut.first);
        return true;
      };
      break;
    case Args::DisplayMode::opengl:
      if (!gl_prepare()) {
        return -1;
      }
      display_func = [](int i, DonutFrame const& donut) {
        glDisplayPoints(donut.second);
        return !glfwWindowShouldClose(g_window);
      };
      break;
  }

  while (true) {
    for (int i{0}; i < its; i++) {
      const auto start_time = std::chrono::high_resolution_clock::now();

      const float inc_sum = inc * i;
      if (not display_func(i, compute_donut(inc_sum, i))) {
        goto terminate;
      }

      sleep_for_frequency_Hz(20, start_time);
    }
  }

terminate:
  glfwTerminate();
  return 0;
}
