#!/bin/bash

docker build -f Dockerfile.ub -t ubuntu1804:base .
docker build -f Dockerfile.r -t ubuntu:r35 .
