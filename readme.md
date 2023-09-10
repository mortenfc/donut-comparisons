# A Donut Revolution Computation

Source:
https://www.a1k0n.net/2011/07/20/donut-math.html

Implemented with varying lightsource (Lx, Ly),
computed in compile time with C++20 constexpr functions,
visualized with OpenGL, set --display=opengl/ascii and
--compute=runtime/compiletime
to compare graphical output and running time between the 4 options.

### Ascii vs OpenGL
![An ASCII Donut Revolution](media/ascii_donut.gif) ![An OpenGL Donut Revolution](media/opengl_donut.gif)

### Compute at startup vs compute while running
![Compute at startup](media/opengl_donut_startuptime.gif) ![Compute while running](media/opengl_donut_runtime.gif)

### Build

```bash
mkdir -p build && cd build; cmake .. && make
```

### Run
./donut --help

### TODO
- Add raytracing and a mirror to opengl rendering
- Optimize render_frame to reduce amount of operations used to be able to compute more frames during compilation
- Improve compile speed by using Bazel
- (Maybe) Make resolution an argument