DROP FUNCTION plv8.d3_slopecontours(gtiff bytea,tresholds int);
CREATE OR REPLACE FUNCTION plv8.d3_slopecontours(gtiff bytea,tresholds int)
RETURNS SETOF JSONB
immutable language plv8
as $$
	const flatten = arr => arr.reduce(
	  (acc, val) => acc.concat(
		Array.isArray(val) ? flatten(val) : val
	  ),
	  []
	);

	var startT = new Date();
	var bytes8 = gtiff;
	var bytes16 = new Uint16Array(bytes8.buffer);
	var tiff = GeoTIFF.parse(bytes16.buffer);
	var image = tiff.getImage(),
	  m = image.getHeight(),
	  n = image.getWidth();
	var rasters = image.readRasters();
	var tiepoint = image.getTiePoints()[0];
	var pixelScale = image.getFileDirectory().ModelPixelScale;
	var geoTransform = [tiepoint.x, pixelScale[0], 0, tiepoint.y, 0, -1*pixelScale[1]];
	var invGeoTransform = [-geoTransform[0]/geoTransform[1], 1/geoTransform[1],0,-geoTransform[3]/geoTransform[5],0,1/geoTransform[5]];
	
	var altData = new Array(image.getHeight());
	for (var j = 0; j<image.getHeight(); j++){
	  altData[j] = new Array(image.getWidth());
	  for (var i = 0; i<image.getWidth(); i++){
		  altData[j][i] = rasters[0][i + j*image.getWidth()];
	  }
	}
	var shadedData = new Array(image.getHeight());
	for (var j = 0; j<image.getHeight(); j++){
	  shadedData[j] = new Array(image.getWidth());
	  for (var i = 0; i<image.getWidth(); i++){
		var gradX, gradY;
		if(i==0) gradX = altData[j][i+1] - altData[j][i];
		else if(i==image.getWidth()-1) gradX = altData[j][i] - altData[j][i-1];
		else gradX = (altData[j][i+1] - altData[j][i])/2 + (altData[j][i] - altData[j][i-1])/2;
	
		if(j==0) gradY = altData[j+1][i] - altData[j][i];
		else if(j==image.getHeight()-1) gradY = altData[j][i] - altData[j-1][i];
		else gradY = (altData[j+1][i] - altData[j][i])/2 + (altData[j][i] - altData[j-1][i])/2;
	
		var slope = Math.PI/2 - Math.atan(Math.sqrt(gradX*gradX + gradY*gradY));
		var aspect = Math.atan2(-gradY, gradX);
	
		shadedData[j][i] = slope;
		/*shadedData[j][i] = Math.sin(altituderad) * Math.sin(slope)
		  + Math.cos(altituderad) * Math.cos(slope)
		  * Math.cos(azimuthrad - aspect);
	    */
	  }
	}
	//plv8.elog(NOTICE,shadedData);
	var contours = d3.contours()
	    .size([n, m])
	    .thresholds(tresholds)
	    .smooth(true)
	    (flatten(shadedData));
	var endT = new Date();
	plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);

	return contours;
$$;

/*EXAMPLE USES:*/
select plv8.plv8_startup();
do language plv8 'load_module("geotiff")';
SET postgis.enable_outdb_rasters = True;
SET postgis.gdal_enabled_drivers = 'ENABLE_ALL';

WITH bounds AS (
	SELECT (ST_MakeEnvelope(120344,488936, 120370,488957,28992)) geom
)
,slope AS (
	SELECT 
	plv8.d3_slopecontours(ST_AsTiff(ST_Clip(ST_Union(rast), geom)),10) AS contour 
	FROM ahn3_raster.rast50cm_tiled, bounds
	WHERE ST_Intersects(rast, geom)
	GROUP BY geom
)
SELECT ST_GeomFromGeoJson(contour::TEXT) FROM slope;
/**/