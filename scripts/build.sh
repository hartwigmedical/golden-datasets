#!/bin/bash

type -a docker > /dev/null

tag_name="eucancan_benchmark"

echo "Building ${tag_name}:latest"
docker build -t "${tag_name}:latest" .