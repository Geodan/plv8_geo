CREATE OR REPLACE FUNCTION public.d3_contour(arr JSON, tresholds int)
RETURNS SETOF JSONB
immutable language plv8
as $$
	const flatten = arr => arr.reduce(
	  (acc, val) => acc.concat(
	    Array.isArray(val) ? flatten(val) : val
	  ),
	  []
	);

	plv8.elog(NOTICE,'NumValues: ' + arr.length);
	var startT = new Date();
	let n = arr.length;
	let m = arr[0].length;
	var contours = d3.contours()
	    .size([n, m])
	    .thresholds(tresholds)
	    .smooth(true)
	    (flatten(arr));
	var endT = new Date();
	plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);
	return contours;
$$;

CREATE OR REPLACE FUNCTION public.d3_contour(arr JSON, tresholds int[])
RETURNS SETOF JSONB
immutable language plv8
as $$
	const flatten = arr => arr.reduce(
	  (acc, val) => acc.concat(
	    Array.isArray(val) ? flatten(val) : val
	  ),
	  []
	);

	plv8.elog(NOTICE,'NumValues: ' + arr.length);
	var startT = new Date();
	let n = arr.length;
	let m = arr[0].length;

	var contours = d3.contours()
	    .size([n, m])
	    .thresholds(tresholds)
	    .smooth(true)
	    (flatten(arr));
	var endT = new Date();
	plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);
	return contours;
$$;



CREATE OR REPLACE FUNCTION public.d3_contour(arr JSON)
  RETURNS JSONB AS
' SELECT public.d3_contour($1,10) '
  LANGUAGE sql STABLE STRICT
  COST 100;





/*EXAMPLE USES:
select plv8_startup();
do language plv8 'load_module("d3")';
do language plv8 'load_module("d3_contour")';


WITH foo AS (
	SELECT ST_SetValue(ST_AddBand(ST_MakeEmptyRaster(3, 3, 0, 0, 1, -1, 0, 0, 0), 1, '8BUI', 1, 0), 1, 2, 5) AS rast
) 
SELECT d3_contour(array_to_json(ST_DumpValues(rast, 1))) AS values FROM foo;

DROP TABLE IF EXISTS tmp.d3contour;
CREATE TABLE tmp.d3contour AS

WITH bounds AS (
	SELECT ST_MakeEnvelope(93784,463805, 93784 + 200,463805 + 200,28992) geom
)


,foo AS (
	SELECT ST_Clip(ST_Union(rast), geom) rast
	FROM ahn3_raster.rasters, bounds
	WHERE ST_Intersects(rast, geom)
	GROUP BY geom
)

,dump AS (
	SELECT array_to_json(ST_DumpValues(rast,1)) AS values, rast 
	FROM foo
	
)
,contours AS (
	SELECT 
		d3_contour(values) contour, rast
	FROM dump
)
SELECT 
	contour->>'value' as z, 
	ST_Translate(
		St_Scale(
			ST_GeomFromGeoJson(contour::TEXT)
			,ST_ScaleX(rast),ST_ScaleY(rast)
		)
		,ST_UpperleftX(rast),ST_UpperleftY(rast)
	) geom
FROM contours,bounds
*/