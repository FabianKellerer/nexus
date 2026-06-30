#!/bin/bash

cmake -DCMAKE_INSTALL_PREFIX=/lhome/ific/f/fkellere/NEWCalibMC/nexus/ -S . -B /lhome/ific/f/fkellere/NEWCalibMC/nexus/build
cmake --build /lhome/ific/f/fkellere/NEWCalibMC/nexus/build --target install -j 8
