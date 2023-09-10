# A Donut Revolution Computation

Source:
https://www.a1k0n.net/2011/07/20/donut-math.html

Implemented with varying lightsource (Lx, Ly),
computed in compile time with C++20 constexpr functions,
visualized with OpenGL, set --display=opengl/ascii and
--compute=runtime/compiletime
to compare graphical output and running time between the 4 options.

![An ASCII Donut Revolution](media/ascii_donut.gif)

![An OpenGL Donut Revolution](media/opengl_donut.gif)

### Build

```bash
mkdir -p build && cd build; cmake .. && make
```

### Run
./donut --help

### TODO
- Add raytracing and a mirror to opengl rendering
- (Maybe) Make resolution an argument for the runtime calculation
- Improve compile speed by using Bazel