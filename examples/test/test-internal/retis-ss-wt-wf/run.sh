#!/usr/bin/env bash
set -e
make clean

folders=('ss-wf-wt-sh' 'sh-ss-wf-wt_cap' 'sh-ss-wt-wf_ha')
moves=("'ss','wf','wt','sh'" "'sh','ss','wf','wt'"
       "'ss','sh','wt','wf'")

for i_f in {0..2}; do
    mkdir -p "${folders[$i_f]}"
    cd "${folders[$i_f]}"
    cp ../initial.xyz .
    cp ../retis.rst .
    cp ../compare.py .
    sed -i "s/'sh', 'sh', 'sh', 'sh'/${moves[$i_f]}/g" retis.rst
    sed -i "s/xxx/${folders[$i_f]}/g" compare.py
    if [ "${folders[$i_f]}" = 'sh-ss-wf-wt_cap' ]; then
        sed -i '40iinterface_cap = -0.55' retis.rst
    fi
    if [ "${folders[$i_f]}" = 'sh-ss-wt-wf_ha' ]; then
        sed -i 's/high_accept = False/high_accept = True/g' retis.rst
    fi
    pyretisrun -i retis.rst -p
    pyretisanalyse -i retis.rst
    python compare.py

    cd ..
done
make clean
