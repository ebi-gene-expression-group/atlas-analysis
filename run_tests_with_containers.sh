#!/usr/bin/env bash

docker run \
  -v $(pwd)/tests:/usr/local/tests \
  --entrypoint=/usr/local/tests/run-tests.bats \
  quay.io/ebigxa/atlas-analysis-base:0.0.1
