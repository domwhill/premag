save_dir='LT2_lineouts'
if [[ ! -d ../OUTPUT/${save_dir} ]]
then
    echo mkdir ../OUTPUT/${save_dir}
    mkdir -p ../OUTPUT/${save_dir}
fi

python2 plot_dqy.py "tot x" ${save_dir}
python2 plot_dqy.py "tot y" ${save_dir}
python2 plot_dqy.py "RL y" ${save_dir}
python2 plot_dqy.py "vN y" ${save_dir}
python2 plot_dqy.py "SH y" ${save_dir}


