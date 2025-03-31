FROM docker:stable

RUN apk add --update bash

COPY run-opensearch.sh /run-opensearch.sh

ENTRYPOINT ["/run-opensearch.sh"]
