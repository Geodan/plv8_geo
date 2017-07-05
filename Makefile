all:
	make js
	make inserts
	make plv8geo--0.0.1.sql

fileclean:
	rm -rf js

js:
	mkdir js
	cd js && \
	wget https://unpkg.com/topojson@3.0.0 -O topojson.js && \
	wget https://d3js.org/d3.v4.min.js -O d3.min.js && \
	wget https://d3js.org/d3-contour.v1.min.js -O d3.contour.min.js


#
inserts:js
	D3:=$(shell cat js/d3.min.js)
	echo "insert into plv8_modules values ('d3',true,$$js$$ $(D3) $$js$$) ON CONFLICT (modname) DO UPDATE SET code = $$js$$ $(D3) $$js$$ WHERE plv8_modules.modname = 'd3';" > sql/ins_d3.sql
	#echo "DELETE FROM plv8_modules WHERE modname = 'topojson';" > sql/ins_topojson.sql
	#echo "COPY plv8_modules (modname, load_on_start, code) FROM stdin;" >> sql/ins_topojson.sql
	#printf "topojson	t	" >> sql/ins_topojson.sql
	#sed ':begin;$!N;s/\n/\\n/;tbegin' < js/topojson.js | sed 's/\t/\\t/g' >> sql/ins_topojson.sql
	#echo "DELETE FROM plv8_modules WHERE modname = 'd3';" > sql/ins_d3.sql
	#printf "INSERT INTO plv8_modules VALUES (d3,t,'" >> sql/ins_d3.sql
	#printf "d3        t       " >> sql/ins_d3.sql
	#sed ':begin;$!N;s/\n/\\n/;tbegin' < js/d3.min.js | sed 's/\t/\\t/g' >> sql/ins_d3.sql
	#printf %q '$(STR)' >> sql/ins_d3.sql
	#cat js/d3.min.js >> sql/ins_d3.sql
	#printf "');" >> sql/ins_d3.sql


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
	cp dump.sql plv8geo--0.0.1.sql
	#cat headerfile.sql > $@
	#cat $^ >> $@

EXTENSION = plv8geo        # the extensions name
DATA = plv8geo--0.0.1.sql  # script files to install
REGRESS = plv8geo_test     # our test script file (without extension)
# postgres build stuff
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
