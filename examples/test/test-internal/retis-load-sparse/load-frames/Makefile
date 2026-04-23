.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find . ! -path ./results ! -path ./results/\* -name \*.txt -delete
	@find -name \out.rst* -delete
	@find -name \*.restart -delete
	@find -name \*.log -delete
	@rm -rf 0*
	@rm -rf pippo/0*
