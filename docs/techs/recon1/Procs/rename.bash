#!/bin/bash
for file in ../Micrographs_part?/*.mrcs ; do ln -sv $file raw${file//[^0-9]/}.mrcs ; done

