DROP FUNCTION IF EXISTS plv8.jsts_voronoi(jsonb);
CREATE OR REPLACE FUNCTION plv8.jsts_voronoi(points jsonb)
  RETURNS jsonb AS
$BODY$
	var geomFact = new jsts.geom.GeometryFactory();
	var reader = new jsts.io.GeoJSONReader();
	var writer = new jsts.io.GeoJSONWriter();
	
	var startT = new Date();
	var runVoronoi = function (sitesWKT) {
	    var sites = reader.read(sitesWKT);
	    plv8.elog(NOTICE, 'reader: ' + (new Date() - startT) / 1000);
	    var builder = new jsts.triangulate.VoronoiDiagramBuilder();
	    builder.setSites(sites);
	    plv8.elog(NOTICE, 'builder: ' + (new Date() - startT) / 1000);
	    var result = builder.getDiagram(geomFact);
	    //plv8.elog(NOTICE, 'diagram: ' + (new Date() - startT) / 1000);
	    
	    return result;
	  }
	var feat = {
		"type": 'Feature',
		"geometry":points
		};
	
	
	var computedTri = runVoronoi(feat.geometry);
	computedTri.normalize();
	var endT = new Date();
	plv8.elog(NOTICE, 'CalcTime: ' + (endT - startT) / 1000);
	return writer.write(computedTri) ;
$BODY$
  LANGUAGE plv8 IMMUTABLE
  COST 100;
