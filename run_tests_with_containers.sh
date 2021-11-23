#!/usr/bin/env bash

docker build -t test/atlas-analysis tests

docker run \
  -v $(pwd):/usr/local/ \
  --entrypoint=/usr/local/tests/run-tests.bats \
  test/atlas-analysis
