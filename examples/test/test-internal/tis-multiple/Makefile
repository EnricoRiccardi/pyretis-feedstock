.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find -name \*out.rst* -delete
	@find -name \*.restart -delete
	@find -name \*.log* -delete
	@rm -f tis-0*
	@rm -rf 0*
	@rm -rf report
	@find results -name "*.html" ! -name "tis-multiple_report.html" -delete
	@rm -f EXIT
