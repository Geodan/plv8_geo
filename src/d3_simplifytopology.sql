
/**
d3_simlifytopology
**/

DROP FUNCTION IF EXISTS plv8.d3_SimplifyTopology(JSONB, numeric);
CREATE FUNCTION plv8.d3_SimplifyTopology(topology JSONB,factor numeric)
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
