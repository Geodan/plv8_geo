DROP FUNCTION  IF EXISTS d3_arctogeom(JSONB, JSONB);
CREATE FUNCTION d3_arctogeom(arc JSONB, transform JSONB)
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


DROP FUNCTION  d3_topobbox(JSONB, JSONB);
CREATE FUNCTION d3_topobbox(arc JSONB, transform JSONB)
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

DROP FUNCTION  d3_ToGeoJson(JSONB,JSONB,JSONB);
CREATE FUNCTION d3_ToGeoJson(entity JSONB,arcs JSONB, transform JSONB)
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


DROP FUNCTION  d3_ToTopojson(collection JSONB, numeric);
CREATE FUNCTION d3_ToTopojson(collection JSONB,q numeric)
RETURNS JSONB
immutable language plv8
as $$
	var startT = new Date();
	var topo = topojson.topology({entities: collection},q);
	var endT = new Date();
	//plv8.elog(NOTICE,'Topotime: ' + (endT - startT)/1000);
	//plv8.elog(NOTICE,JSON.stringify(topo));
	return topo;
$$;

DROP FUNCTION d3_MergeTopology(JSONB, TEXT);
CREATE FUNCTION d3_MergeTopology(topology JSONB,mergekey TEXT)
RETURNS SETOF JSONB
immutable language plv8
as $$
	var startT = new Date();
	var feats = [];
	var data = topojson.feature(topology, topology.objects.entities).features;
	var fips = d3.map(data, function(d){return d.properties[mergekey];}).keys()
	//plv8.elog(NOTICE, JSON.stringify(fips));
	// and merge by fips
	fips.forEach(function(fip) {
	    //var geometries = [topology.objects.entities.geometries[0]];
	    
	    var geometries = topology.objects.entities.geometries.filter(
		function(d) {
			//plv8.elog(NOTICE, JSON.stringify(d));
			return d.properties[mergekey]== fip; 
		}
	    );
	    var geom = topojson.merge(topology, geometries);
	    var feat = {
		type: "Feature", 
		geometry: geom,
		properties: {}
	    };
	    feat.properties[mergekey] = fip;
	    //plv8.elog(NOTICE, JSON.stringify(feat));
	    feats.push(feat);
	});
	var endT = new Date();
	//plv8.elog(NOTICE,'Mergetime: ' + (endT - startT)/1000);
	return feats;
$$;
DROP FUNCTION d3_TopologyToFeatures(topology JSONB);
CREATE FUNCTION d3_TopologyToFeatures(topology JSONB)
RETURNS SETOF JSONB
immutable language plv8
as $$
	var startT = new Date();
	//plv8.elog(NOTICE,JSON.stringify(topology));
	var collection = topojson.feature(topology, topology.objects.entities) 
	var endT = new Date();
	//plv8.elog(NOTICE,'Conversiontime: ' + (endT - startT)/1000);
	return collection.features;
$$;

DROP FUNCTION d3_SimplifyTopology(JSONB, numeric);
CREATE FUNCTION d3_SimplifyTopology(topology JSONB,factor numeric)
RETURNS SETOF JSONB
immutable language plv8
as $$
	var startT = new Date();
	//plv8.elog(NOTICE,JSON.stringify(topology));
	var presimplified = topojson.presimplify(topology);
	var collection = topojson.simplify(presimplified,factor);
	var endT = new Date();
	//plv8.elog(NOTICE,'Conversiontime: ' + (endT - startT)/1000);
	return collection;
$$;