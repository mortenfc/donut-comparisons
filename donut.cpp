#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <array>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>

#include <chrono>
#include <cmath>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <thread>

#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

//* Why do big resolutions + big angle spacing causes black dots to appear?
//* When projecting the 3d torus to 2d as XY, the rotational points won't overlap if it does not compute a point frequently enough,
//* leading to ooz can select a point on another surface between points in the surface facing the viewer.

// WARNING: Changing these rendering constants can heavily affect compilation duration and make it operations overflow during compulation, especially if
// DonutStorage is constexpr. -fconstexpr-ops-limit=173554432

// One-frame donut rendering
constexpr uint32_t screen_width = 500;
constexpr uint32_t screen_height = screen_width;
constexpr uint8_t opengl_resolution_mutliplier = 1080 / screen_height;
constexpr float theta_spacing = 0.07 * (50.0 / screen_width);
constexpr float phi_spacing = 0.02 * (50.0 / screen_width);

// Rendering multiple frames (rotating the donut)
constexpr float frame_frequency = 15.0;
constexpr int rev_count = 4;
constexpr float kA_ls = 0.5;
constexpr float kA = 1.0;
constexpr float kB = 1.5;
constexpr int its = 100;
constexpr float radians_one_cycle = float(rev_count) * M_PI;
constexpr float inc = radians_one_cycle / float(its);

using Ascii = std::array<std::array<char, screen_height>, screen_width>;
using Lums = std::array<std::array<float, screen_height>, screen_width>;
using DonutFrame = std::pair<Ascii, Lums>;
using DonutStorage = std::array<DonutFrame, its>;
/*constexpr*/ DonutStorage render_a_rotating_donut();
//!* NOTE Remove constexpr here to compile larger resolutions. But it will be slow to start instead.
/*constexpr*/ static DonutStorage asciis_and_lums{render_a_rotating_donut()};

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

struct RenderParams {
  float A;
  float B;
  float Lx;
  float Ly;
};

struct RenderCoordinates {
  float L;
  float ooz;
  int xp;
  int yp;
  int ls_x;
  int ls_y;
};

// TODO make this constexpr.. But all vectors should probably be changed to array of preknown size steps_to_thread_end - begin_thread_at_step
/*constexpr*/ std::vector<RenderCoordinates> compute_render_coordinates(RenderParams const& p, int thread_idx = 0,
                                                                        int thread_count = 1) {
  const float cosA = cos(p.A), sinA = sin(p.A);
  const float cosB = cos(p.B), sinB = sin(p.B);
  const float theta_range = 2.0F * M_PI;
  const float phi_range = theta_range;

  const int num_theta_steps = static_cast<int>(std::ceil(theta_range / theta_spacing));
  const int num_phi_steps = static_cast<int>(std::ceil(phi_range / phi_spacing));
  const int total_steps = num_theta_steps * num_phi_steps;

  int steps_per_thread = total_steps / thread_count;
  int steps_to_thread_end = (thread_idx + 1) * steps_per_thread;
  int begin_thread_at_step = thread_idx * steps_per_thread;

  std::vector<RenderCoordinates> results(steps_to_thread_end - begin_thread_at_step);

  // Flattened nested loop
  for (int i = begin_thread_at_step; i < steps_to_thread_end; ++i) {
    int theta_i = i / num_phi_steps;  // 1 1 1 2 2 2 3 3 3
    int phi_i = i % num_phi_steps;    // 1 2 3 1 2 3 1 2 3

    float theta = theta_i * theta_spacing;
    float phi = phi_i * phi_spacing;

    // TODO This math can be heavily optimized
    float costheta = cos(theta), sintheta = sin(theta);
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
    const int xp = (int)(screen_width / 2.0F + K1 * ooz * x);
    const int yp = (int)(screen_height / 2.0F - K1 * ooz * y);

    const int ls_x = (int)(((float)screen_width / 2.0F) + p.Lx * ((float)screen_width / 2.01F));
    const int ls_y = (int)(((float)screen_height / 2.0F) - p.Ly * ((float)screen_height / 2.01F));

    const float L = Lz * (sinA * sintheta + cosA * costheta * sinphi) +
                    p.Ly * (costheta * cosphi * sinB + cosB * (cosA * sintheta - costheta * sinA * sinphi)) +
                    p.Lx * (cosB * costheta * cosphi - sinB * (cosA * sintheta - costheta * sinA * sinphi));
    // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
    // is pointing away from us, so we won't bother trying to plot it.

    results[i - begin_thread_at_step] = {L, ooz, xp, yp, ls_x, ls_y};
  }
  return results;
}

// TODO make this constexpr.. But all vectors should probably be changed to array of preknown size steps_to_thread_end - begin_thread_at_step
/*constexpr*/ DonutFrame map_out_results(std::vector<RenderCoordinates> const& results) {
  auto ascii = ascii_init_copy;
  auto zbuffer = zbuffer_init_copy;
  auto lum = lum_init_copy;

  for (const auto& result : results) {
    // Use result.L, result.ooz, result.xp, result.yp to calculate ascii and lum
    if (result.L > -1) {
      // test against the z-buffer. larger 1/z means the pixel is
      // closer to the viewer than what's already plotted.
      if (result.ooz > zbuffer[result.xp][result.yp]) {
        zbuffer[result.xp][result.yp] = result.ooz;
        int luminance_index = (result.L + 1.F) * 2.2;  // max index is 5
        // luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
        // now we lookup the character corresponding to the
        // luminance and plot it in our ascii:
        ascii[result.xp][result.yp] = ".-=?%@"[luminance_index];
        ascii[result.ls_x][result.ls_y] = 'X';

        lum[result.xp][result.yp] = result.L;
        lum[result.ls_x][result.ls_y] = -10;
      }
    }
  }

  return std::make_pair(ascii, lum);
};

DonutFrame render_frame(RenderParams const& p, uint8_t thread_count) {
  std::vector<std::future<std::vector<RenderCoordinates>>> futures;
  for (uint8_t i = 0; i < thread_count; ++i) {
    // Launch a new thread
    futures.push_back(std::async(std::launch::async, compute_render_coordinates, p, i, thread_count));
  }
  // Join the results
  std::vector<RenderCoordinates> results;
  for (auto& future : futures) {
    auto result = future.get();
    results.insert(results.end(), result.begin(), result.end());
  }

  return map_out_results(results);
}

/*constexpr*/ DonutFrame render_frame(RenderParams const& p) {
  return map_out_results(compute_render_coordinates(p));
}

/*constexpr*/ DonutStorage render_a_rotating_donut() {
  DonutStorage out;

  float A{0}, A_ls{0}, B{0};
  float i{0};
  for (float inc_cum{0.0}; inc_cum < radians_one_cycle; inc_cum += inc) {
    A_ls += inc * kA_ls;
    A += inc * kA;
    B += inc * kB;
    float Lx = cos(A_ls);
    float Ly = sin(A_ls);
    out[i] = render_frame(RenderParams{A, B, Lx, Ly});
    i++;
  }

  return out;
}

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

template <typename Type>
class RunningStatistics {
 public:
  RunningStatistics(Type initial_val) : count{1}, mean{initial_val}, min{mean}, max{mean}, variance{0}, std_dev{0} {}

  void update(Type new_val) {
    const int count_m1 = count;
    count++;
    min = std::min<Type>(min, new_val);
    max = std::max<Type>(max, new_val);
    const int old_mean_diff = (new_val - mean);
    if (count_m1 > 0) {
      variance = (count_m1 / count) * variance + (1.0F / count_m1) * old_mean_diff * old_mean_diff;
      std_dev = sqrt<Type>(variance);
    }
    mean = (mean * count_m1 + new_val) / count;
  }

  void print(std::string name = "RunningStatistics") {
    std::cout << std::endl << name << std::endl;
    std::cout << "Mean: " << mean << "\nMin: " << min << "\nMax: " << max << "\nStdDev: " << std_dev << std::endl;
  }

 private:
  Type count;
  Type mean;
  Type min;
  Type max;
  Type variance;
  Type std_dev;
};

void sleep_for_frequency_Hz(float frequency, std::chrono::time_point<std::chrono::high_resolution_clock> start_time) {

  auto target_cycle_time = std::chrono::duration<double>(1.0 / frequency);
  auto processing_time = std::chrono::high_resolution_clock::now() - start_time;
  auto sleep_time_micros = std::chrono::duration_cast<std::chrono::microseconds>(target_cycle_time - processing_time);

  std::cout << "\nFrequency: " << frequency << "\nProcessing time in micros: "
            << std::chrono::duration_cast<std::chrono::microseconds>(processing_time).count()
            << "\nSleep time in micros: "
            << std::chrono::duration_cast<std::chrono::microseconds>(sleep_time_micros).count() << std::endl;

  static RunningStatistics stats{std::chrono::duration_cast<std::chrono::microseconds>(processing_time).count()};
  stats.update(std::chrono::duration_cast<std::chrono::microseconds>(processing_time).count());
  stats.print("#####  Statistics for processing time [ms]  #####");

  // Sleep for the computed time
  if (sleep_time_micros.count() > 0) {
    std::this_thread::sleep_for(sleep_time_micros);
  }
  std::cout << "\033[2J\033[1;1H";
}

struct Args {
  enum class ProcessingMode { startup_time, run_time, run_time_multi };
  ProcessingMode processing_mode = ProcessingMode::startup_time;

  enum class DisplayMode { ascii, opengl };
  DisplayMode display_mode = DisplayMode::opengl;
};

void print_options() {
  std::cout << "Options:\n"
            << "--display=ascii: Set display mode to ascii\n"
            << "--display=opengl: Set display mode to opengl (default)\n"
            << "--compute=run_time: Set compute mode to run time\n"
            << "--compute=startup_time: Set compute mode to compile time (default)\n"
            << "--compute=run_time_multi: Set compute mode to compile time\n";
}

Args parse_arguments(int argc, char* const argv[]) {
  Args args;
  int opt;
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_options();
      exit(0);
    } else if (strncmp(argv[i], "--display=", 10) == 0) {
      if (strcmp(argv[i] + 10, "ascii") == 0) {
        args.display_mode = Args::DisplayMode::ascii;
      } else if (strcmp(argv[i] + 10, "opengl") == 0) {
        args.display_mode = Args::DisplayMode::opengl;
      } else {
        std::cout << "No such option \"" << argv[i] + 10 << "\"\n";
        print_options();
        exit(1);
      }
    } else if (strncmp(argv[i], "--compute=", 10) == 0) {
      if (strcmp(argv[i] + 10, "run_time") == 0) {
        args.processing_mode = Args::ProcessingMode::run_time;
      } else if (strcmp(argv[i] + 10, "startup_time") == 0) {
        args.processing_mode = Args::ProcessingMode::startup_time;
      } else if (strcmp(argv[i] + 10, "run_time_multi") == 0) {
        args.processing_mode = Args::ProcessingMode::run_time_multi;
      } else {
        std::cout << "No such option \"" << argv[i] + 10 << "\"\n";
        print_options();
        exit(1);
      }
    }
  }

  return args;
}

struct ThreadData {
  std::condition_variable cv_main;
  std::condition_variable cv_worker;
  std::mutex mtx;
  std::atomic<bool> finished{false};
  RenderParams data_in[2];
  DonutFrame data_out[2];
};

void run_time_multi(ThreadData& td) {
  int index{0};

  {
    std::unique_lock<std::mutex> temp_lock{td.mtx};
    td.cv_worker.wait(temp_lock);  // Sync first cycle
  }

  while (!td.finished) {
    // Render frame is so slow it should not wait
    index = (index + 1) % 2;  // Move to the next index in the ring buffer
    td.data_out[index] = render_frame(td.data_in[index], std::thread::hardware_concurrency());  // 1 0 1 0
    td.cv_main.notify_one();  // Notify renderer about new data
  }
}

int main(int argc, char* const argv[]) {
  Args args = parse_arguments(argc, argv);

  ThreadData td;
  std::thread worker;

  DonutFrame runtime_donut_frame;
  std::function<bool(int, DonutFrame const&)> display_func;
  std::function<DonutFrame const&(float, int)> compute_donut;

  switch (args.processing_mode) {
    case Args::ProcessingMode::run_time:
      compute_donut = [&runtime_donut_frame](float inc_sum, int i) -> DonutFrame const& {
        float A_ls = inc_sum * kA_ls;
        float A = inc_sum * kA;
        float B = inc_sum * kB;
        float Lx = cos(A_ls);
        float Ly = sin(A_ls);
        runtime_donut_frame = render_frame(RenderParams{A, B, Lx, Ly});
        return runtime_donut_frame;
      };
      break;
    case Args::ProcessingMode::run_time_multi:
      worker = std::thread(run_time_multi, std::ref(td));
      compute_donut = [&runtime_donut_frame, &td, &worker](float inc_sum, int i) -> DonutFrame const& {
        int index = i % 2;

        float A_ls = inc_sum * kA_ls;
        float A = inc_sum * kA;
        float B = inc_sum * kB;
        float Lx = cos(A_ls);
        float Ly = sin(A_ls);

        {
          std::unique_lock<std::mutex> temp_lock{td.mtx};
          td.data_in[index] = RenderParams{A, B, Lx, Ly};  // 0 1 0 1
          if (i == 0) {
            td.cv_worker.notify_one();  // Sync first cycle
          }
          td.cv_main.wait(temp_lock);  // Sync to the slower thread since there is bidirectional data dependency
        }
        runtime_donut_frame = td.data_out[index];  // 0 1 0 1

        return runtime_donut_frame;
      };
      break;
    case Args::ProcessingMode::startup_time:
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

      sleep_for_frequency_Hz(frame_frequency, start_time);
    }
  }

terminate:
  glfwTerminate();
  if (args.processing_mode == Args::ProcessingMode::run_time_multi) {
    td.finished = true;
    td.cv_main.notify_one();
    td.cv_worker.notify_one();
    worker.join();
  }
  return 0;
}
