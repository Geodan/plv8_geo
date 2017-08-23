FROM geodan/docker-postgres-plv8-postgis:9.6-2.1.0-2.3

RUN buildDependencies="build-essential \
    ca-certificates \
    curl \
    git-core \
    postgresql \
    python-psycopg2 \
    libpq-dev \
    postgresql-server-dev-all" \
  && apt-get update \
  && apt-get install -y --no-install-recommends ${buildDependencies} \
  && git clone https://github.com/Geodan/plv8_geo.git \
  && cd plv8_geo \
  && make \
  && make install \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean \
  && apt-get remove -y  ${buildDependencies} \
  && apt-get autoremove -y \
  && mkdir -p /docker-entrypoint-initdb.d

  COPY ./initdb.sh /docker-entrypoint-initdb.d/postgis.sh
