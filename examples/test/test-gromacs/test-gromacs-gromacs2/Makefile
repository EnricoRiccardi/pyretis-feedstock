.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find -name \*.txt -delete
	@find -name out.rst\* -delete
	@find -name \*.restart -delete
	@find -name rnd.state -delete
	@find -name \*.log* -delete
	@find -name mdout.mdp -delete
	@find -name EXIT -delete
	@find -name \*.tpr -delete
	@rm -rf run-gromacs*/gromacs_inp*
	@rm -rf run-gromacs*/retis-run.rst
	@rm -rf run-gromacs*/0*
	@rm -rf run-gromacs*/out*
	@rm -rf gromacs_input/pyretis.mdp
