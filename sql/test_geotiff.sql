do language plv8 'load_module("geotiff")';

do language plv8 $$
var json_result = plv8.execute( "set bytea_output='escape';SELECT ST_AsTiff(ST_AddBand(ST_MakeEmptyRaster(10, 10, 0, 0, 1, -1, 0, 0, 0),ARRAY[ROW(NULL, '16BUI', 255, 0)]::addbandarg[])) AS rast");
var bytes8 = json_result[0].rast;
var bytes16 = new Uint16Array(bytes8.buffer)
plv8.elog(NOTICE, bytes16);
//var buffer = new ArrayBuffer(bytes16.length);
//bytes16.map(function(i, value){buffer[i] = value});
var tiff = GeoTIFF.parse(bytes16.buffer);
var image = tiff.getImage(),
      values = image.readRasters()[0],
      m = image.getHeight(),
      n = image.getWidth();
plv8.elog(NOTICE, m,n);      
$$;