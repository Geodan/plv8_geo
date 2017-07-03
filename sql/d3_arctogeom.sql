DROP FUNCTION  IF EXISTS plv8.d3_arctogeom(JSONB, JSONB);
CREATE FUNCTION plv8.d3_arctogeom(arc JSONB, transform JSONB)
RETURNS geometry
immutable language plv8
AS $$
 
	function decodeArc(transform, arc) {
	  var x = 0, y = 0;
	 
	  return arc.map(function(position) {
	    position = position.slice();
	    position[0] = (x += position[0]) * transform.scale[0] + transform.translate[0];
	    position[1] = (y += position[1]) * transform.scale[1] + transform.translate[1];
	    return position;
	  });
	}

	var startT = new Date();
	var geojson = {
		type: 'LineString',
		coordinates: decodeArc(transform,arc)
	}
	
	var plan = plv8.prepare( 'SELECT ST_GeomFromGeoJSON($1) as bbox', ['text'] );
	var rows = plan.execute( JSON.stringify(geojson) );
	var endT = new Date();
	//plv8.elog(NOTICE,'Topotime: ' + (endT - startT)/1000);
	return rows[0].bbox;
$$;