#!/bin/bash

source config-0.2-$1
cp config-0.2-$1 config

export OMP_NUM_THREADS=1

./fit_eastern.R >& log.fit_eastern_0.2-$1

rsync -av /tmp/PLScomposition_raw_eastern_0.2-$1.nc paciorek@smeagol:/var/tmp/paleon/composition-eastern-moretaxa/tmp/

rsync -av /tmp/PLScomposition_eastern_0.2-$1_full.nc paciorek@smeagol:/var/tmp/paleon/composition-eastern-moretaxa/output/

rsync -av /tmp/sigma2_eastern_0.2-$1.Rda paciorek@smeagol:/var/tmp/paleon/composition-eastern-moretaxa/output/
