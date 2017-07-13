# plv8_geo
PLV8 functions for geospatial data



## Installation
### Prerequisites
 plv8 [https://github.com/plv8/plv8] 

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


