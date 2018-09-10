#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t synerclust/synerclust:${VERSION}
docker build -t synerclust/synerclust:latest


