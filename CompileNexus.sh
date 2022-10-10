#!/bin/bash

cmake -DCMAKE_INSTALL_PREFIX=/data/nexus/ -S . -B /data/nexus/build
cmake --build build --target install -j 4
