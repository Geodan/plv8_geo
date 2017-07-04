CREATE SCHEMA plv8;
CREATE EXTENSION postgis;
CREATE EXTENSION plv8;
DO $$
  plv8.elog(WARNING, 'plv8.version = ' + plv8.version); // Will output the PL/v8 installed as a PostgreSQL `WARNING`.
$$ LANGUAGE plv8;
CREATE EXTENSION plv8geo;
select plv8.plv8_startup();
--do language plv8 'load_module("delaunator")';
--SELECT plv8.delaunator([[1,1],[10,10],[5,5]]) AS values FROM foo;
