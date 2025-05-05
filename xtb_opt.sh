#!/bin/bash

# 脚本使用说明
usage() {
    echo "Usage: $0 <input_file> <data_dir> <start_line> <end_line> <out_dir> <ncores> [split]"
    echo "Parameters:"
    echo "  <input_file>    输入文件路径"
    echo "  <data_dir>      数据文件夹路径"
    echo "  <start_line>    起始行号"
    echo "  <end_line>      终止行号（小于等于0时默认处理从起始行到文件末尾）"
    echo "  <out_dir>       输出文件夹路径"
    echo "  <ncores>        核心数（可选，默认为6）"
    echo "  [split]         子文件夹划分大小（可选，默认为-1，表示不划分）"
    exit 1
}

# 检查参数数量是否正确
if [ $# -lt 5 ]; then
    usage
fi

# 1. 接收参数
input_file=$1
data_dir=$2
start_line=$3
end_line=$4
out_dir=$5
ncores=${6:-6}  # 如果未提供ncores，默认为6
split=${7:--1}   # 如果未提供split，默认为1（表示不划分）

# 转化路径为绝对路径
input_file=$(realpath "$input_file")
data_dir=$(realpath "$data_dir")
out_dir=$(realpath "$out_dir")
xtb_path=$(realpath "$(dirname "$0")/xtb")

# 2. 创建输出文件夹
mkdir -p "$out_dir" || { echo "无法创建输出文件夹 $out_dir"; exit 1; }

# 3. 设置环境变量
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS="$ncores,1"
export MKL_NUM_THREADS="$ncores"

# 4. 创建临时文件夹
temp_dir=$(mktemp -d temp_XXXXXXXXXXXXXXXX)
cd "$temp_dir" || { echo "无法切换到临时目录 $temp_dir"; exit 1; }

# 计算总行数
total_lines=$(wc -l < "$input_file")

# 处理起始行和终止行
start_line=$((start_line < 1 ? 1 : start_line))
if [ "$end_line" -le 0 ]; then
    end_line=$total_lines
else
    end_line=$((end_line > total_lines ? total_lines : end_line))
fi

realline=$(($end_line - $start_line + 1))

# 如果split>0，创建所有可能的子文件夹
if [ "$split" -gt 0 ]; then
    # 预先创建所有可能的子文件夹
    for ((i=0; i<total_lines; i+=split)); do
        start=$((i + 1))
        end=$((i + split))
        # if [ "$end" -gt "$total_lines" ]; then
        #     end=$total_lines
        # fi
        mid_path="${start}-${end}"
        mkdir -p "$out_dir/$mid_path" || { echo "无法创建子文件夹 $out_dir/$mid_path"; exit 1; }
    done
fi

# 记录开始时间
start_time=$(date +%s)

# 处理各行
for ((current_line=start_line; current_line<=end_line; current_line++)); do
    # 获取当前行内容
    line=$(sed -n "${current_line}p" "$input_file")
    
    # 解析前三个部分
    a=$(echo "$line" | cut -f1)
    b=$(echo "$line" | cut -f2)
    c=$(echo "$line" | cut -f3)

    # echo ${a} ${b} ${c}
    # ask user whether go on
    # exit

    # 提取文件编号
    num=$(echo "$c" | grep -oE '[0-9]+' || echo "0")
    
    # 计算mid_path
    if [ "$split" -gt 0 ]; then
        tdx=$(( (num - 1) / split ))  # 整除
        start=$((tdx * split + 1))
        end=$((tdx * split + split))
        mid_path="${start}-${end}"
    else
        mid_path=""
    fi
    
    # 复制文件
    if [ "$split" -gt 0 ]; then
        cp "$data_dir/$mid_path/$c" "conf.xyz" || { echo "无法复制 $data_dir/$mid_path/$c"; continue; }
    else
        cp "$data_dir/$c" "conf.xyz" || { echo "无法复制 $data_dir/$c"; continue; }
    fi
    
    # 运行xtb命令
    ${xtb_path} conf.xyz -c "$b" --uhf 0  -o > tempfile.txt 2>&1
    cont=$(cat tempfile.txt)
    
    # 检查是否收敛
    if grep -q "Property Printout" tempfile.txt; then
        if [ "$split" -gt 0 ]; then
            cp "xtbopt.xyz" "$out_dir/$mid_path/$c" || { echo "无法复制 $temp_dir/xtbopt.xyz"; }
        else
            cp "xtbopt.xyz" "$out_dir/$c" || { echo "无法复制 $temp_dir/xtbopt.xyz"; }
        fi
    fi
    
    # 输出进度
    if [ $(( (current_line - start_line + 1) % 1000 )) -eq 0 ]; then
        elapsed=$(( $(date +%s) - $start_time  ))
        echo "=============== $(($current_line - $start_line + 1)) / $realline done , $elapsed seconds used ============="
    fi
done

# 输出最终进度
elapsed=$(( $(date +%s) - $start_time  ))
echo "=============== $realline / $realline done , $elapsed seconds used ============="

# 删除临时文件夹
cd - || { echo "无法切换回原目录"; }
rm -rf "$temp_dir"
