# A Donut Comparison

Source:
https://www.a1k0n.net/2011/07/20/donut-math.html

Implemented with varying lightsource (Lx, Ly),
computed in compile time with C++20 constexpr functions,
visualized with OpenGL. 
Compare startup computation, single threaded run-time and slight multithreading.
C++20 can't compute these huge data structures at compile-time on my computer because it takes too much memory apparently.

### Display comparison: Ascii vs OpenGL
![An ASCII Donut Revolution](media/ascii_donut.gif) ![An OpenGL Donut Revolution](media/opengl_donut.gif)

### Compute comparisons
![Compute at startup](media/opengl_donut_startup_time_stats.gif) 

![Compute while running single threaded](media/opengl_donut_run_time_stats.gif)

![Compute while running slightly multi threaded](media/opengl_donut_run_time_multi_stats.gif)

// TODO compute while running the rendering multithreaded.

### Results for OpenGL for splitting rendering and visualization into 1 thread each
| Processing time [ms] 	| Mean    	|
|----------------------	|---------	|
| run_time             	| ~110000 	|
| run_time_multi       	| ~100000 	|
| startup_time         	| ~10000  	|

Doing the rendering in a separate thread is only slightly faster, since visualizing with opengl only takes about ~1/10 of the time.
The rendering itself should be split into multiple threads.

### Comparison for rendering in all available threads
Lower resolution required after rewriting program to use less heap or it memory overflows.
Should just go back to commit 917eafd11a1e159b802f127a2eb6d00e5f4585eb and reserve

| Mean processing time [ms] | Debug    	| Release  |
|--------------------------	|---------	| ---------|
| run_time                	| ~20057 	|          |
| run_time_multi        	| ~14391 	|          |
| startup_time          	| ~6099  	|          |

### Build

```bash
mkdir -p build && cd build; cmake .. && make -j
```

### Run
./donut --help

### TODO
- Update results with rendering being split into many threads
- Add raytracing and a mirror to opengl rendering
- Optimize render_frame to reduce amount of operations used to be able to compute more frames during compilation
- Improve compile speed by using Bazel
- (Maybe) Make resolution an argumen