#!/bin/bash

GCM_TEST_PATH=..
GCM_TEST=$GCM_TEST_PATH/crystal_plasticity_grains

if [ $# -lt 2 ]; then
  echo "usage: run_intgration_alg.sh [NP] [tag]"
  echo "exit"
else
  NP=$1
  tag=$2

  np=1

  out_dir=out_$tag
  rm -rf $out_dir
  mkdir $out_dir
  cd $out_dir

  mpirun -np $NP $GCM_TEST $GCM_TEST_PATH/params_$tag.in $GCM_TEST_PATH/CO/co $np
fi
