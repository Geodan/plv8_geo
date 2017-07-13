select plv8.plv8_startup();
do language plv8 ' load_module("d3"); ';
do language plv8 ' load_module("d3-geo"); ';
do language plv8 ' load_module("d3-force"); ';
do language plv8 ' load_module("d3-hexbin"); ';
do language plv8 ' load_module("d3-contour"); ';
do language plv8 ' load_module("delaunator"); ';
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