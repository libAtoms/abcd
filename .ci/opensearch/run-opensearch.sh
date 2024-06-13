#!/usr/bin/env bash
source ./.ci/opensearch/functions/imports.sh
set -euxo pipefail

if [[ -z $OPENSEARCH_VERSION ]]; then
  echo -e "\033[31;1mERROR:\033[0m Required environment variable [OPENSEARCH_VERSION] not set\033[0m"
  exit 1
fi

OPENSEARCH_REQUIRED_VERSION="latest"
# Starting in 2.12.0, security demo configuration script requires an initial admin password
if [ "$OPENSEARCH_VERSION" != "$OPENSEARCH_REQUIRED_VERSION" ]; then
  OPENSEARCH_INITIAL_ADMIN_PASSWORD="admin"
fi

for (( node=1; node<=${NODES-1}; node++ ))
do
  port=$((PORT + $node - 1))

  if [[ "$SECURITY_ENABLED" == "true" ]]; then
    healthcmd="curl -vvv -s --insecure -u admin:$OPENSEARCH_INITIAL_ADMIN_PASSWORD --fail https://localhost:$port/_cluster/health || exit 1"
    security=($(cat <<-END

END
  ))
  elif [[ "$SECURITY_ENABLED" == "false" ]]; then
    healthcmd="curl -vvv -s --fail http://localhost:$port/_cluster/health || exit 1"
    security=($(cat <<-END
      --env plugins.security.disabled=true
END
  ))
  fi

  docker run \
    --rm \
    --detach \
    --name="os${node}" \
    --env "cluster.name=docker-opensearch" \
    --env "http.port=${port}" \
    --env discovery.type=single-node \
    --env bootstrap.memory_lock=true \
    --env "OPENSEARCH_JAVA_OPTS=-Xms4g -Xmx4g" \
    --env OPENSEARCH_INITIAL_ADMIN_PASSWORD=$OPENSEARCH_INITIAL_ADMIN_PASSWORD \
    "${security[@]}" \
    --publish "${port}:${port}" \
    --ulimit nofile=65536:65536 \
    --ulimit memlock=-1:-1 \
    --health-cmd="$(echo $healthcmd)" \
    --health-interval=2s \
    --health-retries=20 \
    --health-timeout=2s \
    opensearchproject/opensearch:${OPENSEARCH_VERSION}

  if wait_for_container "os$node"; then
    echo -e "\033[32;1mSUCCESS:\033[0m OpenSearch up and running\033[0m"
  fi
done
