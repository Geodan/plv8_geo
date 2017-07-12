
all:
	make src/inserts.sql
	make plv8geo--0.0.1.sql

cleanfiles:
	rm src/inserts.sql
	rm plv8geo--0.0.1.sql
	
src/inserts.sql:
	printf "/** This will insert all javascript modules into a table **/" > src/inserts.sql
	printf "insert into plv8_modules values ('topojson',false,\$$js\$$" >> src/inserts.sql 
	cat js/topojson.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql
	
	printf "insert into plv8_modules values ('delaunator',false,\$$js\$$" >> src/inserts.sql 
	cat js/delaunator.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql
	
	printf "insert into plv8_modules values ('d3',false,\$$js\$$" >> src/inserts.sql
	cat js/d3.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql
	
	printf "insert into plv8_modules values ('d3-contour',false,\$$js\$$" >> src/inserts.sql
	cat js/d3-contour.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql
	
	printf "insert into plv8_modules values ('d3-force',false,\$$js\$$" >> src/inserts.sql
	cat js/d3-force.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql
	
	printf "insert into plv8_modules values ('d3-hexbin',false,\$$js\$$" >> src/inserts.sql
	cat js/d3-hexbin.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql
	
	printf "insert into plv8_modules values ('d3-geo',false,\$$js\$$" >> src/inserts.sql
	cat js/d3-geo.js >> src/inserts.sql
	printf  "%s\$$js\$$);" >> src/inserts.sql

functions := src/inserts.sql \
	src/d3_arctogeom.sql \
	src/d3_contour.sql \
	src/d3_hexbin.sql \
	src/d3_mergetopology.sql \
	src/d3_simplifytopology.sql \
	src/d3_slopecontours.sql \
	src/d3_togeojson.sql \
	src/d3_topobbox.sql \
	src/d3_topologytofeatures.sql \
	src/d3_totopojson.sql \
	src/delaunator.sql \
	src/slope.sql 

plv8geo--0.0.1.sql: $(functions) 
	cat src/headerfile.sql > $@
	cat $^ >> $@

EXTENSION = plv8geo        # the extensions name
DATA = plv8geo--0.0.1.sql  # script files to install
REGRESS = plv8geo_test     # our test script file (without extension)
# postgres build stuff
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
