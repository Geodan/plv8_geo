
/**
d3_contour
**/
DROP FUNCTION IF EXISTS plv8.d3_contour(JSON, int);
CREATE OR REPLACE FUNCTION plv8.d3_contour(gtiff bytea, tresholds int)
RETURNS SETOF JSONB
immutable language plv8
as $$
	var startT = new Date();
	var bytes8 = gtiff;
	var bytes16 = new Uint16Array(bytes8.buffer);
	//FIXME: namespace for geotiff
	var tiff = GeoTIFF.parse(bytes16.buffer);
	var image = tiff.getImage(),
      values = image.readRasters()[0],
      m = image.getHeight(),
      n = image.getWidth();
	plv8.elog(NOTICE,'Width/Height ', n,m);
	plv8.elog(NOTICE,'NumValues: ' + values.length);
	
	var contours = d3.contours()
	    .size([n, m])
	    .thresholds(tresholds)
	    .smooth(true)
	    (values);
	var endT = new Date();
	plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);
	return contours;
$$;

DROP FUNCTION IF EXISTS plv8.d3_contour(gtiff bytea, tresholds int[]);
CREATE OR REPLACE FUNCTION plv8.d3_contour(gtiff bytea, tresholds int[])
RETURNS SETOF JSONB
immutable language plv8
as $$

	var startT = new Date();
	var bytes8 = gtiff;
	var bytes16 = new Uint16Array(bytes8.buffer);
	//FIXME: namespace for geotiff
	var tiff = GeoTIFF.GeoTIFF.parse(bytes16.buffer);
	var image = tiff.getImage(),
      values = image.readRasters()[0],
      m = image.getHeight(),
      n = image.getWidth();
	plv8.elog(NOTICE,'Width/Height ', n,m);
	plv8.elog(NOTICE,'NumValues: ' + values.length);

	var contours = d3.contours()
	    .size([n, m])
	    .thresholds(tresholds)
	    .smooth(true)
	    (values);
	var endT = new Date();
	plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);
	return contours;
$$;

/*TODO fix this
DROP FUNCTION IF EXISTS plv8.d3_contour(arr JSON);
CREATE OR REPLACE FUNCTION plv8.d3_contour(arr JSON)
  RETURNS JSONB AS
' SELECT plv8.d3_contour($1,10) '
  LANGUAGE sql STABLE STRICT
  COST 100;*/





/*
EXAMPLE USES:
select plv8.plv8_startup();
do language plv8 'load_module("d3")';
do language plv8 'load_module("d3-contour")';


WITH foo AS (
	SELECT ST_SetValue(ST_AddBand(ST_MakeEmptyRaster(3, 3, 0, 0, 1, -1, 0, 0, 0), 1, '8BUI', 1, 0), 1, 2, 5) AS rast
) 
SELECT d3_contour(array_to_json(ST_DumpValues(rast, 1))) AS values FROM foo;




do language plv8 'load_module("d3")';
do language plv8 'load_module("d3-contour")';
do language plv8 'load_module("geotiff")';

DROP TABLE IF EXISTS tmp.tmp;
CREATE TABLE tmp.tmp AS
WITH bounds AS (
	SELECT ST_MakeEnvelope(137236,467884 , 143010,472643,28992) geom
)
,args AS (
	SELECT ROW(1, '-10-300:-10-300', '8BUI', NULL)::reclassarg arg
)
,contours AS (
	SELECT plv8.d3_contour(ST_AsTiff(
		ST_Reclass(
		  ST_Resample(
			ST_Clip(ST_Union(rast), geom)
		  ,10,10)
		,arg)
		),ARRAY[0,5,10,15,20,25,30]) contour ,
	ST_Clip(ST_Union(rast),geom) as rast
	FROM ahn3_raster.rasters, bounds, args
	WHERE ST_Intersects(rast, geom)
	GROUP BY geom,arg
)
SELECT 
	(contour->>'value')::double precision as z, 
	ST_Translate(
		St_Scale(
			ST_GeomFromGeoJson(contour::TEXT)
			,ST_ScaleX(rast),ST_ScaleY(rast)
		)
		,ST_UpperleftX(rast),ST_UpperleftY(rast)
	) geom
FROM contours

*/