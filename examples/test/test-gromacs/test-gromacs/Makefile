.PHONY: clean

clean:
	@find -name \*.pyc -delete
	@find -name \*.pyo -delete
	@find -name __pycache__ -delete
	@find -name \*.txt -delete
	@find -name out.rst\* -delete
	@find -name \*.restart -delete
	@find -name rnd.state -delete
	@find -name \*.log -delete
	@find -name mdout.mdp -delete
	@find -name pyretis.mdp -delete
	@find -name EXIT -delete
	@find -name \*.tpr -delete
	@rm -rf 0-forward-step
	@rm -rf 0-backward-step
	@rm -rf 0-backward-plain
	@rm -rf 0-forward-plain
	@rm -rf 1-forward-step
	@rm -rf 1-backward-step
	@rm -rf 1-backward-plain
	@rm -rf 1-forward-plain
	@rm -rf 2-forward-step
	@rm -rf 2-backward-step
	@rm -rf 2-backward-plain
	@rm -rf 2-forward-plain
	@rm -f engine-run.rst
