
/**
d3_hexbin
**/

CREATE OR REPLACE FUNCTION plv8.d3_hexbin(arr1 JSON, arr2 JSON,radius integer)
RETURNS SETOF JSONB
immutable language plv8
as $$
if (arr1.length != arr2.length){
	plv8.elog(ERROR, 'Arrays are not same length');
}
let arr = arr1.map((d,i)=>{
	return {0:d[0],1:d[1],"key":arr2[i]};
});

//plv8.elog(NOTICE,JSON.stringify(arr));
var startT = new Date();
var hexbin = d3.hexbin()
    .extent([[d3.min(arr1, d=>d.x), d3.min(arr1, d=>d.y)], [d3.max(arr1, d=>d.x), d3.max(arr1, d=>d.y)]])
    .radius(radius);
let bins = hexbin(arr);
let res = bins.map(b=>{
	//Silly map because plv8 has problems with mixed arrays
	return {x:b.x,y:b.y, data: b};
});
//plv8.elog(NOTICE,JSON.stringify(res));
var endT = new Date();
plv8.elog(NOTICE,'CalcTime: ' + (endT - startT)/1000);

return res;

$$;
/* TEST 
select plv8.plv8_startup();
do language plv8 'load_module("d3")';
do language plv8 'load_module("d3-hexbin")';

SELECT plv8.d3_hexbin(('[[1,2],[0.5,0.5],[2,2]]')::json,'["foo","bar","baz"]'::JSON,1);

DROP TABLE IF EXISTS tmp.hexbin;
CREATE TABLE tmp.hexbin AS 
WITH 
bounds AS (
SELECT ST_Transform(ST_MakeEnvelope(524596,6855645 ,562126,6876847,3857),28992) box
)
,palen AS (
	SELECT 
	ARRAY[ST_X(geom), ST_Y(geom)] geom
	,'test'::text as val
	,gebrksdoel 
	FROM bagagn_201704.adressen , bounds WHERE ST_Intersects(box,geom)
--AND gebrksdoel = 'onderwijsfunctie'
)
,hexbin AS (
	SELECT plv8.d3_hexbin(array_to_json(array_agg(geom)),array_to_json(array_agg(val)),200) AS value 
	,gebrksdoel 
	FROM palen
	GROUP BY gebrksdoel 
)

SELECT ST_MakePoint((value->>'x')::double precision, (value->>'y')::double precision,28992) geom,
jsonb_array_length(value->'data') num 
,gebrksdoel 
FROM hexbin;

*/
