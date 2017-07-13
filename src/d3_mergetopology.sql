
/**
d3_mergetopology
**/

DROP FUNCTION IF EXISTS plv8.d3_MergeTopology(JSONB, TEXT);
CREATE FUNCTION plv8.d3_MergeTopology(topology JSONB,mergekey TEXT)
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
