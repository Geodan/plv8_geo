
/**
d3_topobbox
**/
DROP FUNCTION IF EXISTS plv8.d3_topobbox(JSONB, JSONB);
CREATE FUNCTION plv8.d3_topobbox(arc JSONB, transform JSONB)
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
	var xmin = d3.min(arc,d=>d[0]);
	var xmax = d3.max(arc,d=>d[0]);
	var ymin = d3.min(arc,d=>d[1]);
	var ymax = d3.max(arc,d=>d[1]);
	/* TODO: future work to get a binary geometry back
	var wkx = require('wkx');
	var wkt = 'POLYGON(('+xmin +' '+ymin+','+xmax+' '+ ymax+'))';
	var geometry = wkx.Geometry.parse(wkt);
	*/
	var plan = plv8.prepare( 'SELECT ST_MakeEnvelope($1,$2,$3,$4) as bbox', ['float','float','float','float'] );
	var rows = plan.execute( xmin,ymin, xmax, ymax );
	var endT = new Date();
	//plv8.elog(NOTICE,'Topotime: ' + (endT - startT)/1000);
	return rows[0].bbox;
$$;
