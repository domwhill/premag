#!/bin/bash
environment=$1
envexists=$(conda env list | grep ${environment})
if [ "${envexists}" == "" ]; then
	echo "no environment ${environment} found"
	conda env create -f environment.yml
else
	echo "environmant ${environment} already exists. Nothing to be done."
fi
