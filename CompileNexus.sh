#!/bin/bash

sudo cmake -DCMAKE_INSTALL_PREFIX=/Users/fabian/Documents/nexus/ -S . -B /Users/fabian/Documents/nexus/build
sudo cmake --build build --target install -j 8
