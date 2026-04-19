.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find -name \*.txt -delete
	@find -name \out.rst* -delete
	@find -name \*.restart -delete
	@find -name \*.log* -delete
	@find -name compare.png -delete
	@find -name report -type d -exec rm -rf {} + 2>/dev/null; true
