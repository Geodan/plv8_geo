CREATE OR REPLACE FUNCTION plv8.slope(gtiff bytea)
RETURNS JSONB
immutable language plv8
as $$
	var startT = new Date();
	var bytes8 = gtiff;
	var bytes16 = new Uint16Array(bytes8.buffer);
	var tiff = GeoTIFF.parse(bytes16.buffer);
	var image = tiff.getImage();
	
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
	var endT = new Date();
	plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);
	//plv8.elog(NOTICE,shadedData);
	return shadedData;
$$;

/*EXAMPLE USES:
select plv8.plv8_startup();
do language plv8 'load_module("d3")';
do language plv8 'load_module("d3_contour")';


WITH foo AS (
	SELECT ST_SetValue(ST_AddBand(ST_MakeEmptyRaster(3, 3, 0, 0, 1, -1, 0, 0, 0), 1, '8BUI', 1, 0), 1, 2, 5) AS rast
) 
SELECT plv8.slope(ST_AsTiff(rast)) AS values FROM foo;

*/