DROP FUNCTION plv8.d3_TopologyToFeatures(topology JSONB);
CREATE FUNCTION plv8.d3_TopologyToFeatures(topology JSONB)
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