.PHONY:
	clean

clean:
	find -name \*.pyc -delete
	find -name \*.pyo -delete
	find -name __pycache__ -delete
	find -name EXIT -delete
	find -name out.rst\* -delete
	find -name \*.restart -delete
	find -name rnd.state -delete
	find -name \*.log -delete
	find -name mdout.mdp -delete
	find -name \*.tpr -delete
	-rm -rf results/pyvisa_compressed_data.hdf5
	-rm -rf pyvisa_compressed_data*
	-rm -rf 0*
	-rm -rf *_0*
	-rm -rf retis-load-rc-run.rst
	-rm -rf pippo/0*
	-rm -rf gromacs_input/pyretis.mdp
    
