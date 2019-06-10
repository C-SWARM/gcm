#!/usr/bin/env bash

exe=$1
ref=$2
out=$3
compare="numdiff --absolute-tolerance=1.0e-14 --relative-tolerance=1.0e-14 $ref $out"

# cleanup on failure.
function cleanup {
    echo "Removing $out"
    rm -f $out
}
trap cleanup EXIT

echo $exe
$exe
echo $compare
$compare
