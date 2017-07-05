\echo Use "CREATE EXTENSION plv8geo" to load this file. \quit
create table IF NOT EXISTS plv8_modules(modname text primary key, load_on_start boolean, code text);
CREATE SCHEMA IF NOT EXISTS plv8geo;

