all:
	make inserts
	make plv8geo--0.0.1.sql

inserts:
	echo 'test' 

	
	#delaunator:=$(shell cat js/delaunator.js)
	#echo "insert into plv8_modules values ('d3',true,$$js$$ $(delaunator) $$js$$) > sql/ins_delaunator.sql


inserts:js
	printf 'insert into plv8_modules values (''d3'',true,$$js$$' > sql/ins_d3.sql
	cat js/d3.min.js >> sql/ins_d3.sql
	printf '$$js$$);' >> sql/ins_d3.sql
	printf 'insert into plv8_modules values (''topojson'',true,$$js$$' > sql/ins_topojson.sql
	cat js/topojson.js >> sql/ins_topojson.sql
	printf '$$js$$);' >> sql/ins_topojson.sql


functions := sql/ins_topojson.sql \
	sql/ins_d3.sql \
	sql/d3_arctogeom.sql \
	sql/d3_contour.sql \
	sql/d3_hexbin.sql \
	sql/d3_mergetopology.sql \
	sql/d3_simplifytopology.sql \
	sql/d3_slopecontours.sql \
	sql/d3_togeojson.sql \
	sql/d3_topobbox.sql \
	sql/d3_topologytofeatures.sql \
	sql/d3_totopojson.sql \
	sql/delaunator.sql \
	sql/slope.sql 

plv8geo--0.0.1.sql: $(functions) 
	#cp dump.sql plv8geo--0.0.1.sql
	cat headerfile.sql > $@
	cat $^ >> $@

EXTENSION = plv8geo        # the extensions name
DATA = plv8geo--0.0.1.sql  # script files to install
REGRESS = plv8geo_test     # our test script file (without extension)
# postgres build stuff
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
