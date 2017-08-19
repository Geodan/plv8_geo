/**
earcut
**/
DROP FUNCTION IF EXISTS plv8.earcut(geom JSONB);
CREATE OR REPLACE FUNCTION plv8.earcut(geom JSONB)
RETURNS JSONB
immutable language plv8
as $$
var startT = new Date();
var data = earcut.flatten(geom.coordinates);
var triangles = earcut(data.vertices); 
var endT = new Date();
//plv8.elog(NOTICE, 'CalcTime: ' + (endT - startT) / 1000);
return triangles;
$$;

/* TEST:
SELECT plv8.earcut(ST_AsGeoJson(ST_MakeEnvelope(0,0,10,10))::JSONB);
*/