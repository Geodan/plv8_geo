/**
delaunator
**/
DROP FUNCTION IF EXISTS plv8.delaunator(arr numeric[]);
CREATE OR REPLACE FUNCTION plv8.delaunator(points JSONB)
RETURNS JSONB
immutable language plv8
as $$
var points = points.coordinates;
var startT = new Date();
var delaunay = new Delaunator(points);
var endT = new Date();
//plv8.elog(NOTICE, 'CalcTime: ' + (endT - startT) / 1000);
return delaunay.triangles;
$$;


/** AGGREGATOR FUNCTION

Comment: Unforunately, the aggregated version is a LOT slower than the direct triangulation with GeoJSON
**/

DROP TYPE IF EXISTS plv8.dpoint CASCADE;
CREATE TYPE plv8.dpoint AS(x numeric, y numeric, z numeric);

DROP FUNCTION IF EXISTS plv8.delaunator_state(points numeric[], point numeric[]) CASCADE;
CREATE OR REPLACE FUNCTION plv8.delaunator_state(points plv8.dpoint[], point plv8.dpoint)
RETURNS plv8.dpoint[] AS
$$
//plv8.elog(NOTICE,'state',points);
if (points == null) points = [];
if (point == null) return [points];
else {
	points.push(point);
	//plv8.elog(NOTICE,'state returning ',JSON.stringify(points));
	return points;
}
$$ language plv8 immutable;

--DROP FUNCTION IF EXISTS plv8.delaunator_final(points numeric[]) CASCADE;
CREATE OR REPLACE FUNCTION plv8.delaunator_final(points plv8.dpoint[])
RETURNS JSONB AS
$$
//plv8.elog(NOTICE,'final',points);
var coords = points.map(d => [d.x, d.y]);
var coordinates = [];
var delaunay = new Delaunator(coords);
let triangles = delaunay.triangles;
for (var i = 0; i < triangles.length; i += 3) {
	coordinates.push([
		points[triangles[i]],
		points[triangles[i + 1]],
		points[triangles[i + 2]]
	]);
}
return coordinates;

$$ language plv8 immutable;

--DROP AGGREGATE IF EXISTS plv8.delaunator_agg(numeric[]);
CREATE AGGREGATE plv8.delaunator_agg(plv8.dpoint)(
	SFUNC = plv8.delaunator_state,
	STYPE = plv8.dpoint[], FINALFUNC = plv8.delaunator_final

);
/*TODO:
working on a way to get normals for every triangle
*/
CREATE OR REPLACE FUNCTION plv8.delaunator_normals(triangles JSONB)
RETURNS JSONB AS
$$

function triangleNormal(p0, p1, p2, output) {
	if (!output) output = []
	let x0 = p0.x;
	let y0 = p0.y;
	let z0 = p0.z;

	let x1 = p1.x;
	let y1 = p1.y;
	let z1 = p1.z;

	let x2 = p2.x;
	let y2 = p2.y;
	let z2 = p2.z;

	var p1x = x1 - x0
	var p1y = y1 - y0
	var p1z = z1 - z0

	var p2x = x2 - x0
	var p2y = y2 - y0
	var p2z = z2 - z0

	var p3x = p1y * p2z - p1z * p2y
	var p3y = p1z * p2x - p1x * p2z
	var p3z = p1x * p2y - p1y * p2x

	var mag = Math.sqrt(p3x * p3x + p3y * p3y + p3z * p3z)
	if (mag === 0) {
		output[0] = 0
		output[1] = 0
		output[2] = 0
	} else {
		output[0] = p3x / mag
		output[1] = p3y / mag
		output[2] = p3z / mag
	}

	return output
}
let normals = triangles.map(t => {
	return triangleNormal(t[0], t[1], t[2]);
});
//plv8.elog(NOTICE, 'normals for', JSON.stringify(normals));

$$ language plv8 immutable;

/*EXAMPLE USES:
select plv8.plv8_startup();
do language plv8 'load_module("delaunator")';
WITH points AS ( 
	SELECT PC_Explode(pa) pt FROM ahn3_pointcloud.vw_ahn3 
	WHERE PC_Intersects(ST_MakeEnvelope(120339,488721, 120749,489032,28992),pa)
	LIMIT 1000000
)-- 1.1 seconds (1,000,000 points) SELECT count(*) FROM points
,delaunator_agg AS (
	SELECT 
	plv8.delaunator_agg((PC_Get(pt,'x'),PC_Get(pt,'y'))::plv8.dpoint)
	FROM points
)
,delaunator AS (
	SELECT plv8.delaunator(ST_AsGeoJson(ST_Collect(Geometry(pt)))::JSONB)
	FROM points
) --9.6 seconds (1.556 calctime)
,stdelaunay AS (
	SELECT ST_DelaunayTriangles(ST_Collect(Geometry(pt))) triags
	FROM points
) --17.6 seconds
,triangulate2dz AS (
	SELECT ST_Triangulate2dz(ST_Collect(Geometry(pt))) triags
	FROM points
) --13.1 seconds
SELECT count(*) FROM delaunator




select plv8.plv8_startup();
do language plv8 'load_module("delaunator")';
SELECT plv8.delaunator([[1,1],[10,10],[5,5]]) AS values FROM foo;
*/