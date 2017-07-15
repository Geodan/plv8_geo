# plv8_geo
A postgis extension with PLV8 functions for geospatial data

This extension will load a couple of hand-picked javascript geometry-related libraries into your postgres database to use them with plv8. At the moment the following libraries are included (also see the js directory):

d3: "^4.7.4",
d3-contour: "^1.1.0",
d3-force: "^1.0.6",
d3-geo: "^1.6.3",
d3-hexbin: "^0.2.2",
delaunator: "git+https://github.com/tomvantilburg/delaunator.git",
geotiff: "^0.4.1",
topojson": "^3.0.0"

## Installation
### Prerequisites
 PostgreSQL 6.4 with PL/v8 2.0.3 [https://github.com/plv8/plv8] 
 note: PL/v8 2.0.3 packages not available yet, requires manual build See the instructions over here: https://github.com/plv8/plv8/blob/master/doc/plv8.md#installing-plv8 (make sure to use make static in order to get the latest v8 engine)

### Build
On a linux machine (tested on Ubuntu 17.04), clone the git repository and run:
```
make
sudo make install
```

## Testing
You can test to see if the library will load and run on your system with:
```
make installcheck (optional: PGUSER=username PGHOST=myhost etc..)
```
### Creating extension
In your sql prompt run `CREATE EXTENSION plv8geo`
This will put all the plv8geo stuff into a new schema called plv8. Functions will be available from this schema.

## List of functions

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


