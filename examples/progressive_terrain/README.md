# Progressive Terrain Mesh Example

![For example](../../resources/screenshots/terrain/hex_terrain6.png?raw=true)

Download/extract/convert shuttle radar topography mission data with 1 arc-second resolution ("SRTM1") to visualize with a progressive mesh.

Data is publicly available per [USGS](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects).

Use `prepare_srtm1.py` to generate PNG formatted grayscale height maps (~5.5 GB for all PNGs specified by default options).
Use `min/max_lat/lon` parameters to download only a subset.

Use `view_terrain.py` to visualize the PNG height maps using a progressive mesh which can be subdivided for better resolution.
Use the `LEFT` and `RIGHT` arrow keys to increase/decrease the resolution.
Use the `UP` and `DOWN` arrow keys to increase/decrease the sealevel.
Use `F2` to toggle rendering of edges (or `F1` to toggle rendering of faces).

![An example using ~24k triangles](screenshots/24576-tris.png?raw=true "24k triangles")

Using a modest laptop from 2014 (Dell Latitude E7440), ~400k triangle meshes can be slow to load into memory but are managable to view.
