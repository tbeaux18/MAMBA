#!/bin/bash

docker build -f Dockerfile.ubuntubase -t ubuntu1804:base .
docker build -f Dockerfile.rbase -t ubuntu:r35 .
docker build -f Dockerfile.fastqc -t ubuntu:fastqc .
