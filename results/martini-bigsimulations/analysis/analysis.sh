#!/bin/bash

mkdir rdf
mkdir clusters
mkdir crystals
mkdir tilt
mkdir neighbors

for DIR in $( ls .. | grep 'f11-SNda' ); do

    sed "s/VARDIR/$DIR/g" martini_cnt_cluster_analysis/rdf.sge       > rdf-$DIR.sge
    sed "s/VARDIR/$DIR/g" martini_cnt_cluster_analysis/clusters.sge  > clusters-$DIR.sge
    sed "s/VARDIR/$DIR/g" martini_cnt_cluster_analysis/crystals.sge  > crystals-$DIR.sge
    sed "s/VARDIR/$DIR/g" martini_cnt_cluster_analysis/tilt.sge      > tilt-$DIR.sge
    sed "s/VARDIR/$DIR/g" martini_cnt_cluster_analysis/neighbors.sge > neighbors-$DIR.sge

    qsub rdf-$DIR.sge
    qsub clusters-$DIR.sge
    qsub crystals-$DIR.sge
    qsub tilt-$DIR.sge
    qsub neighbors-$DIR.sge

done

