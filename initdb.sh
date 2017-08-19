#!/bin/sh

set -e

# Perform all actions as $POSTGRES_USER
export PGUSER="$POSTGRES_USER"

"${psql[@]}" <<- 'EOSQL'
    ALTER USER postgres WITH PASSWORD 'postgres';
    \connect postgres
    CREATE EXTENSION IF NOT EXISTS plv8geo cascade;
\q
EOF