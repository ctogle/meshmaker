# Progressive Terrain Mesh Example

Download/extract/convert shuttle radar topography mission data with 1 arc-second resolution ("SRTM1") to visualize with a progressive mesh.
Data is publicly available per [USGS](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects).

![Example](../../resources/screenshots/terrain/hex_terrain6.png?raw=true)

Use `prepare_srtm1.py` to generate PNG formatted grayscale height maps (~5.5 GB for all PNGs specified by default options).

 - Use `min/max_lat/lon` parameters to download only a subset.

Use `view_terrain.py` to visualize the PNG height maps using a progressive mesh which can be subdivided for better resolution.

 - Use `LEFT` and `RIGHT` arrow keys to increase/decrease resolution.
 - Use `UP` and `DOWN` arrow keys to increase/decrease sealevel.
 - Use `F2` to toggle edge rendering (or `F1` to toggle face rendering).

![Example](screenshots/6-tris.png?raw=true "6 triangles")

![Example](screenshots/96-tris.png?raw=true "~96 triangles")

![Example](screenshots/1536-tris.png?raw=true "~1.5k triangles")

![Example](screenshots/24576-tris.png?raw=true "~24k triangles")

Using a modest laptop from 2014 (Dell Latitude E7440), a ~390k triangle mesh can be slow to load into memory but is managable to view.

![Example](screenshots/393216-tris.png?raw=true "393k triangles")
