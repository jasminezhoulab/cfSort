#!/bin/bash

likelihood_ratio_threshold=2

mkdir -p output

gzip -dc input/example.reads.txt.gz | ../src/cfsort -a em.global.unknown -p $likelihood_ratio_threshold -O tissueFraction+readCountPerBillion -r stdin -o output/example.class.specific.read.counts.profile -T 7 -t input/example.markers.txt

echo ""
echo "==============="
cmp -s output/example.class.specific.read.counts.profile output/example.class.specific.read.counts.profile_reference && echo "cfsort Test run passed" || echo "cfsort Test run failed"
echo "==============="
