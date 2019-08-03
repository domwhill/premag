#!/bin/bash
in_string="site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')"
out_string="site.addsitedir('/Users/'+ userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')"
#in_string="test"
#out_string="test2"
echo ${in_string}
echo sed -i '' "s#${in_string}#${out_string}#g" **/*.py
sed -i '' "s#${in_string}#${out_string}#g" **/*.py
