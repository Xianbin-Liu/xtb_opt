# intro
`rdkit ETKDGv3`来计算初始化构象以及[xtb](https://xtb-docs.readthedocs.io/en/latest/sp.html)方法来优化分子构象。

# how to use
## installation
1. 预编译版本：https://github.com/grimme-lab/xtb/releases/latest
直接下载解压，并把xtb-dist/bin/添加到环境变量中。
2. python环境
```bash
pip install ase xtb pandas rdkit tqdm
```
## 获取初始构象
```bash
python prepare_conf.py
```
## 优化分子构象
```bash
./xtb_opt.sh <input_file> <data_dir> <start_line> <end_line> <out_dir> <ncores> [split]

./xtb_opt.sh data/amlist_test.txt data/conf_test 0 0 data/conf_test_opt
```
`amlist_test`作为input file，每一行有3个`\t`分割的项：需要优化的smiles、charges、对应的初始conf_path。
`conf_test`为初始conf_path的目录，`conf_test_opt`为优化后的conf_path的目录。
`start_line`和`end_line`为需要优化的分子在input file中的行号，`ncores`为使用的核数，`split`为是否分割成多个文件进行优化。

`xtb_opt.sh`中的核心代码为
```bash
${xtb_path} conf.xyz -c "$b" --uhf 0  -o > tempfile.txt 2>&1
```

# prepare conf
将`.`连接的smiles分割，分别作canonical得到所有需要初始化构象的smiles（去除atom mapping以减少重复），smiles2conf保存cano smiles到conf path的映射关系。

# load conf

`utils.py`

load conf和prepare过程保持一致，对整个`input smiles`按`.`分割 - 取得cano的smiles后计算带am的“有效原子”个数 - 计算不同部分的connection mask - 取得conf - （如果需要addHs） 将xyz以及分部分load的分子，对它们的原子次序rerank成`MolFromSmiles(input smiles)`的次序