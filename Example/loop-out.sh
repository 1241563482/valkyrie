
#! /bin/bash
mkdir -p output
for i in `ls POSCAR-* | sed 's/POSCAR-//g'`
do
cp $i/vasprun.xml output/"vasprun.xml-"$i
done
phonopy -f output/*
phonopy -s -p band.conf --tolerance=0.1
phonopy-bandplot --gnuplot > ph.out
#phonopy -t mesh.conf --tolerance=0.1
#python get_free_energy.py
#phonopy --dim="2 2 2" --irreps="0 0 0 1e-3"
