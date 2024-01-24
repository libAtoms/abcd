#!/usr/bin/env bash
source ./.ci/opensearch/functions/imports.sh
set -euxo pipefail

if [[ -z $OPENSEARCH_VERSION ]]; then
  echo -e "\033[31;1mERROR:\033[0m Required environment variable [OPENSEARCH_VERSION] not set\033[0m"
  exit 1
fi

for (( node=1; node<=${NODES-1}; node++ ))
do
  port=$((PORT + $node - 1))

  if [[ "$SECURITY_ENABLED" == "true" ]]; then
    healthcmd="curl -vvv -s --insecure -u admin:admin --fail https://localhost:$port/_cluster/health || exit 1"
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
