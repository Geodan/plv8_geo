EXTENSION = plv8geo        # the extensions name
DATA = plv8geo--0.0.1.sql  # script files to install
REGRESS = plv8_test     # our test script file (without extension)

# postgres build stuff
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)