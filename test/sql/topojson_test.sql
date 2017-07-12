select plv8.plv8_startup();
do language plv8 ' load_module("topojson"); ';

/**
do language plv8 $$
	var topo = topojson.topology({entities: {
		type: 'Feature',
		geometry: {
			type: 'Point',
			coordinates: [0,0]
		}
	}},1);
	plv8.elog(NOTICE,JSON.stringify(topo));
$$;
**/