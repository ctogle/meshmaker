# meshmaker
Framework for procedural mesh generation/visualization

The ultimate goal of this repository is a pipeline for generating content suitable for games, 
specifically large cities with unique, fully decorated building interiors, realistic roads, terrain, foliage, etc.

There are already suitable [open-source](https://github.com/magnificus/Procedural-Cities) starting points 
for such an endeavor as well as commercial products that cover the same domain ([Houdini](https://www.sidefx.com/)), 
but this repository is for self-education / fun; 
it is written almost entirely in Python without leaning heavily on external packages where possible 
and would be between these two examples in scope if in a completed state.

Currently it contains, to varying degrees of "prototyped", basic geometric primitives used in computational geometry 
(3-vector, quaternion, 2-D fields, planar polygons, planar graphs), 
classes for manipulation/visualization 
(2-manifold mesh, [BSP trees](https://github.com/evanw/csg.js/), 
[delaunay triangulation](https://people.eecs.berkeley.edu/~jrs/meshpapers/delnotes.pdf), L-systems, 
shader based OpenGL pipeline using [ModernGL](https://github.com/moderngl/moderngl), scenegraphs, OBJ-exporter), 
and other supporting methods/classes for generating content using these objects 
(schema for arbitrary parameterization of a 2-manifold surface, 
[medial axis approximation based on triangulation](https://www.sciencedirect.com/science/article/abs/pii/S0965997812000828), 
terrain generation from 
[SRTM](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects) 
data, L-systems).
