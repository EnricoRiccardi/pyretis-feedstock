.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find -name \*.txt -delete
	@find -name \*.wfn -delete
	@find -name \*out.rst* -delete
	@find -name \*.restart -delete
	@find -name \*.log -delete
	@find -name \*.bak-* -delete
	@rm -rf 0*
	@rm -f cp2k_input/pyretis.inp
	@rm -rf backward
	@rm -rf backward-plain
	@rm -rf forward
	@rm -rf forward-plain
	@rm -rf genvel
	@rm -rf forward-single-step
