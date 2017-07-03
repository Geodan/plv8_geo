
/**
d3_togeojson
**/
DROP FUNCTION IF EXISTS plv8.d3_ToGeoJson(JSONB,JSONB,JSONB);
CREATE FUNCTION plv8.d3_ToGeoJson(entity JSONB,arcs JSONB, transform JSONB)
RETURNS JSONB
immutable language plv8
as $$
	var startT = new Date();
	
	var result = [];

	function decodeArc(transform, arc) {
	  var x = 0, y = 0;
	  return arc.map(function(position) {
	    position = position.slice();
	    position[0] = (x += position[0]) * transform.scale[0] + transform.translate[0];
	    position[1] = (y += position[1]) * transform.scale[1] + transform.translate[1];
	    return position;
	  });
	}
	
	var arcsundelta = arcs.map(function(a){
		return {id: a.id, coords: decodeArc(transform, a.data)};
		//return {id: a.id, coords: a.data};
	});
	
	entity.arcs[0].forEach(function(id,i){
		//watch out, negative index is -1 based
		var coords = arcsundelta.find(d=>(d.id == id && id >=0)||(d.id == Math.abs(id)-1 && id <0) ).coords;
		if (id < 0){
		  coords.reverse();
		}
		if (i<entity.arcs[0].length){
			coords.pop();//remove last element
		}
		coords.forEach(c=>result.push(c));
	});
	result.push(result[0]);
	delete entity.arcs;
	entity.geometry = {
		coordinates: [result],
		type: entity.type
	};
	entity.type = 'Feature';
	
	var endT = new Date();
	//plv8.elog(NOTICE,'Topotime: ' + (endT - startT)/1000);
	//plv8.elog(NOTICE,JSON.stringify(entity));	
	return entity;
$$;
