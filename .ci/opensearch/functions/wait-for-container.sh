#!/usr/bin/env bash
#
# Exposes a routine scripts can call to wait for a container if that container set up a health command
#
# Please source .ci/functions/imports.sh as a whole not just this file
#
# Version 1.0.1
# - Initial version after refactor
# - Make sure wait_for_contiainer is silent
# From https://github.com/opensearch-project/opensearch-py/blob/main/.ci/functions/wait-for-container.sh

function container_running {
  if [[ "$(docker ps -q -f name=$1)" ]]; then
    return 0;
    else return 1;
  fi
}

function wait_for_container {
  set +x
  until ! container_running "$1" || (container_running "$1" && [[ "$(docker inspect -f "{{.State.Health.Status}}" ${1})" != "starting" ]]); do
    echo ""
    docker inspect -f "{{range .State.Health.Log}}{{.Output}}{{end}}" ${1}
    echo -e "\033[34;1mINFO:\033[0m waiting for node $1 to be up\033[0m"
    sleep 4;
  done;

  # Always show logs if the container is running, this is very useful both on CI as well as while developing
  if container_running $1; then
    docker logs $1
  fi

  if ! container_running $1 || [[ "$(docker inspect -f "{{.State.Health.Status}}" ${1})" != "healthy" ]]; then
    echo -e "\033[31;1mERROR:\033[0m Failed to start $1 in detached mode beyond health checks\033[0m"
    echo -e "\033[31;1mERROR:\033[0m dumped the docker log before shutting the node down\033[0m"
    return 1
  else
    echo
    echo -e "\033[32;1mSUCCESS:\033[0m Detached and healthy: ${1}\033[0m"
    return 0
  fi
}
