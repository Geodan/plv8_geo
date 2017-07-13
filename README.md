# plv8_geo
PLV8 functions for geospatial data



## Installation
### Prerequisites
 plv8 [https://github.com/plv8/plv8] 
 plv8_libloader [https://github.com/geodan/plv8_libloader]

### Build
On a linux machine (tested on Ubuntu 17.04), clone the git repository and run:
```
make
sudo make install
```
### Creating extension
In your sql prompt run `CREATE EXTENSION plv8geo`
This will put all the plv8geo stuff into a new schema called plv8. Functions will be available from this schema.

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
SELECT plv8.d3_contour(array_to_json(ST_DumpValues(rast, 1))) AS values FROM foo;
```


