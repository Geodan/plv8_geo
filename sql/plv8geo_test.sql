CREATE SCHEMA plv8;
CREATE EXTENSION postgis;
CREATE EXTENSION plv8;
DO $$
  plv8.elog(WARNING, 'plv8.version = ' + plv8.version); // Will output the PL/v8 installed as a PostgreSQL `WARNING`.
$$ LANGUAGE plv8;
CREATE EXTENSION plv8geo;
create or replace function plv8.plv8_startup()
returns void
language plv8
as
$$
load_module = function(modname) {
  var rows = plv8.execute("SELECT code from plv8_modules " +" where modname = $1", [modname]);
  for (var r = 0; r < rows.length; r++) {
	    var code = rows[r].code;
	    eval("(function() { " + code + "})")();
	  }      
	};
	$$;

	select plv8.plv8_startup(); 
select plv8.plv8_startup();
do language plv8 ' load_module("topojson"); ';
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
