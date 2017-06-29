# plv8_geo
PLV8 functions for geospatial data

## Installation
Make sure you first install plv8 and the plv8_libloader [https://github.com/geodan/plv8_libloader]

Just run the individual SQL files when needed.

## List of functions

### d3_contour
Usage: 
d3_contour()
```sql
select plv8_startup();
do language plv8 'load_module("d3")';
do language plv8 'load_module("d3_contour")';


WITH foo AS (
	SELECT ST_SetValue(ST_AddBand(ST_MakeEmptyRaster(3, 3, 0, 0, 1, -1, 0, 0, 0), 1, '8BUI', 1, 0), 1, 2, 5) AS rast
) 
SELECT d3_contour(array_to_json(ST_DumpValues(rast, 1))) AS values FROM foo;
```


