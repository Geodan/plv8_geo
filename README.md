# plv8_geo
A postgis extension with PLV8 functions for geospatial data

This extension will load a couple of hand-picked javascript geometry-related libraries into your postgres database to use them with plv8. At the moment the following libraries are included (also see the js directory):

 - d3: 4.7.4,
 - d3-contour: 1.1.0,
 - d3-force: 1.0.6,
 - d3-geo: 1.6.3,
 - d3-hexbin: 0.2.2,
 - delaunator: 1.0.2,
 - geotiff: 0.4.1,
 - topojson: 3.0.0

## Docker support

Run database with everything installed on Docker:

```
$ docker run -p 5432:5432 geodan/postgis_plv8_geo
```

login with: postgres/postgres

Command for building image:

```
$ docker build -t geodan/postgis_plv8_geo .
```

## Manual Installation
### Prerequisites
 PostgreSQL 9.4 with PL/v8 2.0.3 [https://github.com/plv8/plv8] 
 
 note: PLv8 2.0.3 packages not available yet, requires manual build See the instructions over here: https://github.com/plv8/plv8/blob/master/doc/plv8.md#installing-plv8 (make sure to use make static in order to get the latest v8 engine)
 Older plv8 packages do not support newer ECMAscript version and therefor most libraries do not work (even refuse to install)
 
### Build
On a linux machine (tested on Ubuntu 17.04), clone the git repository and run:
```
make
sudo make install
```
Note: the make step is only needed when you want to add libraries yourself. It would be necessary to edit the Makefile in order for these libraries to be loaded.

### Test
You can test to see if the library will load and run on your system with:
```
make installcheck (optional: PGUSER=username PGHOST=myhost etc..)
```
### Create extension
In your sql prompt run `CREATE EXTENSION plv8geo`
This will put all the plv8geo stuff into a new schema called plv8. Functions will be available from this schema.

In order to have the libraries loaded at startup time, add the following to postgresql.conf: 
`plv8.startproc = 'plv8.plv8_startup'`


## Available functions

### plv8.d3_totopojson
Usage:
```sql
SELECT plv8.d3_totopojson(<featurecollection>::JSONB)
```
Returns topology
### plv8.d3_simplifytopology
Usage:
```sql
SELECT plv8.d3_simplifytopology(<topology>::JSONB,<factor>::numeric)
```
Returns topology::JSONB

### plv8.d3_mergetopology
Usage:
```sql
SELECT plv8.mergetopology(<topojson>::JSONB,<propertykey>::TEXT)
```
Merges the features ins the topology based on a common property.
Return set of geojson::JSONB features

### plv8.d3_topologytofeatures
Usage:
```sql
SELECT plv8.d3_topologytofeatures(<topojson>::JSONB)
```
Returns set of geojson features


### plv8.d3_contour
Usage: 
```sql
WITH foo AS (
	SELECT ST_SetValue(ST_AddBand(ST_MakeEmptyRaster(3, 3, 0, 0, 1, -1, 0, 0, 0), 1, '8BUI', 1, 0), 1, 2, 5) AS rast
) 
SELECT plv8.d3_contour(array_to_json(ST_DumpValues(rast, 1))) AS values FROM foo;
```

### plv8.d3_hexbin
Usage:
```sql
SELECT plv8.d3_hexbin(('[[1,2],[0.5,0.5],[2,2]]')::json,'["foo","bar","baz"]'::JSON,1);
```
returns setof JSONB with {x:x, y:y, data:[]}

### plv8.delaunator
Usage:
```sql
WITH points AS (
	SELECT ST_GeneratePoints(ST_MakeEnvelope(0,0,100,100),10000) geom
)
SELECT plv8.delaunator(ST_AsGeoJson(geom)::JSONB)
FROM points
```
returns JSONB

## Examples

Simplify an existing set of geometries topologically

```sql

```

#### Create contours out of a raster
```sql
select plv8.plv8_startup();
do language plv8 'load_module("d3-contour")';
do language plv8 'load_module("geotiff")';
SET postgis.gdal_enabled_drivers = 'ENABLE_ALL';


WITH foo AS (
	SELECT ST_SetValue(ST_AddBand(ST_MakeEmptyRaster(3, 3, 0, 0, 1, -1, 0, 0, 0), 1, '8BUI', 1, 0), 1, 2, 5) AS rast
) 
,args AS (
	SELECT ROW(1, '-10-300:-10-300', '16BUI', NULL)::reclassarg arg
)
SELECT plv8.d3_contour(
ST_AsTiff(
ST_Reclass(
rast,arg
)
),10) AS values FROM args,foo;
```
