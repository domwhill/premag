ENVIRONMENT=pop_env
RUN_ENV=source activate ${ENVIRONMENT}

python_files=$(shell git ls-files *.py)

all: install_environment build_plots

install_environment:
	cd environment;\
	./install_environment.sh ${ENVIRONMENT}

run_tests:
	$(RUN_ENV);\
	python -m pytest *.py

fix_code_formatting:
	$(RUN_ENV);\
	yapf -i --style=setup.cfg $(python_files)

build_plots:
	./generate_paper_plots.sh
