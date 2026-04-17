.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find -name \*.txt -delete
	@find -name \*.wfn -delete
	@find -name \*out.rst* -delete
	@find -name \*.restart -delete
	@find -name \*log* -delete
	@find -name \*.bak-* -delete
	@rm -rf 0*
	@rm -rf *pyretis-cp2k*
	@rm -rf step.inp
	@rm -rf continue.inp
	@rm -rf conf.xyz
	@rm -f cp2k_input/pyretis.inp
	@rm -rf cp2k-run
