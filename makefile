ENVIRONMENT=pop_env
RUN_ENV=source activate ${ENVIRONMENT}

python_files=$(shell git ls-files *.py)

all: install_environment check_data build_plots

# check that simulation data has been added to correct folder
check_data:
	if [ -d "data" ];\
	then echo "data directory exists";\
	else echo "No data directory found. Please copy simulation data to folder entitled data/ in root of repo" && exit 1;\
	fi

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
