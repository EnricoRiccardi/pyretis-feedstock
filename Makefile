.PHONY:	
	tests
	tests-silent
	tests3
	tests-silent3
	clean
	upload-docs
	coverage
	pydocstyle
	pycodestyle

coverage:
	pytest --cov=pyretis --cov-report=html --cov-report=xml test/
 
tests:
	pytest -v test/

tests-silent:
	pytest test/

tests3:
	pytest -v test/

tests-silent3:
	pytest test/

pydocstyle:
	pydocstyle --count ./pyretis

pycodestyle:
	pycodestyle

clean:
	find -name \*.pyc -delete
	find -name \*.pyo -delete
	find -name __pycache__ -delete
	find -name \*.so -delete

upload-docs:
	scp -r docs/_build/html/* pyretisweb:WWW/
