all:
	make js
	make inserts

fileclean:
	rm -rf js

js:
	mkdir js
	cd js && \
	wget https://unpkg.com/topojson@3.0.0 -O topojson.js 

TOPOJSON := $(shell cat js/topojson.js)
inserts:js
	echo "DELETE FROM plv8_modules WHERE modname = 'topojson';" > sql/ins_topojson.sql
	echo "COPY plv8_modules (modname, load_on_start, code) FROM stdin;" >> sql/ins_topojson.sql
	printf "topojson	t	" >> sql/ins_topojson.sql
	tr -d '\n' < js/topojson.js | tr -d '\t' >> sql/ins_topojson.sql
functions := sql/d3_arctogeom.sql \
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

