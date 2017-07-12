
all:
	make sql/inserts.sql
	make plv8geo--0.0.1.sql

cleanfiles:
	rm sql/inserts.sql
	rm plv8geo--0.0.1.sql
	
sql/inserts.sql:
	printf "/** This will insert all javascript modules into a table **/" > sql/inserts.sql
	printf "insert into plv8_modules values ('topojson',false,\$$js\$$" >> sql/inserts.sql 
	cat js/topojson.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql
	
	printf "insert into plv8_modules values ('delaunator',false,\$$js\$$" >> sql/inserts.sql 
	cat js/delaunator.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql
	
	printf "insert into plv8_modules values ('d3',false,\$$js\$$" >> sql/inserts.sql
	cat js/d3.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql
	
	printf "insert into plv8_modules values ('d3-contour',false,\$$js\$$" >> sql/inserts.sql
	cat js/d3-contour.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql
	
	printf "insert into plv8_modules values ('d3-force',false,\$$js\$$" >> sql/inserts.sql
	cat js/d3-force.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql
	
	printf "insert into plv8_modules values ('d3-hexbin',false,\$$js\$$" >> sql/inserts.sql
	cat js/d3-hexbin.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql
	
	printf "insert into plv8_modules values ('d3-geo',false,\$$js\$$" >> sql/inserts.sql
	cat js/d3-geo.js >> sql/inserts.sql
	printf  "%s\$$js\$$);" >> sql/inserts.sql

functions := sql/inserts.sql \
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
	cat headerfile.sql > $@
	cat $^ >> $@

EXTENSION = plv8geo        # the extensions name
DATA = plv8geo--0.0.1.sql  # script files to install
REGRESS = plv8geo_test     # our test script file (without extension)
# postgres build stuff
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
