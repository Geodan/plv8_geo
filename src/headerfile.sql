\echo Use "CREATE EXTENSION plv8geo" to load this file. \quit
CREATE TABLE IF NOT EXISTS public.plv8_modules(modname text primary key, load_on_start boolean, code text);
CREATE SCHEMA IF NOT EXISTS plv8;

create or replace function plv8.plv8_startup()
returns void
language plv8
as
$$
load_module = function(modname) {
 var rows = plv8.execute("SELECT code from public.plv8_modules " +" where modname = $1", [modname]);
 for (var r = 0; r < rows.length; r++) {
    var code = rows[r].code;
    eval("(function() { " + code + "})")();
  }      
};
$$;

select plv8.plv8_startup(); 
