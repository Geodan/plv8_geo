CREATE EXTENSION plv8geo;
select plv8.plv8_startup();
do language plv8 'load_module("delaunator")';
SELECT plv8.delaunator([[1,1],[10,10],[5,5]]) AS values FROM foo;