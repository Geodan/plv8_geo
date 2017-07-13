select plv8.plv8_startup();
do language plv8 ' load_module("topojson"); ';
WITH entities AS (
	SELECT St_MakeEnvelope(0,0,10,10) geom
	UNION ALL
	SELECT St_MakeEnvelope(0,10,0,20) geom
)
,geometry as (
	SELECT st_asgeojson(geom) geom 
	from entities
)
,features AS (
    select '{"type": "Feature"}'::JSONB ||
       	   jsonb_set('{}'::JSONB,'{geometry}',geom::JSONB)
	as feat
	FROM geometry g 
)
SELECT 
plv8.d3_topologytofeatures(
plv8.d3_totopojson(
	'{"type": "FeatureCollection"}'::JSONB ||
	 jsonb_set('{}'::JSONB, '{features}',jsonb_agg(feat))
	 ,1e8
)
)
FROM features;