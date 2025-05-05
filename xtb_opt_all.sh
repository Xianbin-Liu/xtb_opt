source /data2/xianbin/anaconda3/etc/profile.d/conda.sh
conda activate xtb

data_dir="data"
for part in train val test
do
    input_file="$data_dir/amlist_${part}.txt"
    echo "Processing $input_file"
    # list all files in the dir: ${data_dir}/conf_${part}
    if [ ! -d "${data_dir}/conf_${part}" ]; then
        echo "Directory ${data_dir}/conf_${part} does not exist. Creating it."
        continue
    fi

    echo "starting to process  ${input_file},  ${data_dir}/conf_${part}"
    ./xtb_opt.sh ${input_file} ${data_dir}_conf_${part} 0 0

done