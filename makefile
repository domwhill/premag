ACTIVATE_ENVIRONMENT=pipenv shell
RUN_ENV=source activate pop_env

python_files=$(shell git ls-files *.py)

install_environment:
	cd environment;\
	conda env create -f environment.yml

run_tests:
	$(RUN_ENV);\
	python -m pytest *.py

fix_code_formatting:
	$(RUN_ENV);\
	yapf -i --style=setup.cfg $(python_files)

build_plots:
	./generate_paper_plots.sh