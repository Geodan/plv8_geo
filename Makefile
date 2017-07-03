EXTENSION = plv8geo        # the extensions name

files := sql/d3_arctogeom.sql \
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

plv8geo--0.0.1.sql: $(files) 
	echo '\echo Use "CREATE EXTENSION plv8geo" to load this file. \quit' > $@
	cat $^ >> $@

DATA = plv8geo--0.0.1.sql  # script files to install
REGRESS = plv8_test     # our test script file (without extension)
# postgres build stuff
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

