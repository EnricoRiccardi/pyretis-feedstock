#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
export MPLBACKEND=Agg
cp -r --update=none ../../test-gromacs/test-load/test-load-sparse/load-traj/* .
cp ../../test-gromacs/test-load/test-load-sparse/load-traj/run.sh run_g.sh


sed -i '$ d' run_g.sh
sed -i '$ d' run_g.sh
sed -i '$ d' run_g.sh

echo "sed -i 's/method = load/method = kick/g' retis-load-rc-run.rst" >> run_g.sh
echo "pyretisrun -i retis-load-rc-run.rst -p" >> run_g.sh

./run_g.sh

pyvisa -i retis-load-rc-run.rst -cmp
cp ../lib/compare_*.py .
unzip pyvisa_compressed_data.hdf5.zip 
unzip results/pyvisa_compressed_data.hdf5.zip -d results
python compare_hdf5.py
rm results/pyvisa_compressed_data.hdf5

echo ' ' >> retis-load-rc-run.rst
echo 'Collective-variable' >> retis-load-rc-run.rst
echo '-------------------' >> retis-load-rc-run.rst
echo 'class = Distance' >> retis-load-rc-run.rst
echo 'module = orderp.py' >> retis-load-rc-run.rst
echo 'idx1 = 0' >> retis-load-rc-run.rst
echo 'idx2 = 4' >> retis-load-rc-run.rst

pyvisa -i retis-load-rc-run.rst -recalculate -cmp
find . -mindepth 1 -not -name 'run.sh' -not -name 'results' -not -path './results/pyvisa_compressed_data.hdf5.zip' -delete
