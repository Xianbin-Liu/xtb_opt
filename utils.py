import os
from ase.io import read as aseread
from tqdm import tqdm
from rdkit import Chem
from os.path import join as path_join
import numpy as np

def cano(x):
    return Chem.MolToSmiles(Chem.MolFromSmiles(x))

def clear_map_number(x):
    '''
        clear the atom map number. and cano the smiles
    '''
    mol = Chem.MolFromSmiles(x)
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')
    tmp_x = Chem.MolToSmiles(mol)
    return Chem.MolToSmiles(Chem.MolFromSmiles(tmp_x))

def cano_readd(x):
    '''
    cano the smiles, re-make atom map number [AtomMapNumber=AtomIdx+1]
    '''
    x = clear_map_number(x)
    mol = Chem.MolFromSmiles(x)
    for idx, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(idx + 1)
    return Chem.MolToSmiles(mol)

def fill_all_amap(rxn):
    am_reacs, am_prods = rxn.split('>>')
    am_prods_mol = Chem.MolFromSmiles(am_prods)
    am_reacs_mol = Chem.MolFromSmiles(am_reacs)
    amap_nums = [x.GetAtomMapNum() for x in am_prods_mol.GetAtoms()]
    max_amap = max(amap_nums) if len(amap_nums) > 0 else 0
    for idx, atom in enumerate(am_reacs_mol.GetAtoms()):
        if not atom.hasProp('molAtomMapNumber'):
            max_amap += 1
            atom.SetAtomMapNum(max_amap)
    return Chem.MolToSmiles(am_reacs_mol) + '>>' + am_prods

def load_smiles2conf(data_dir, part):
    """
    Load the smiles to conf mapping from the specified directory and part.

    Args:
        data_dir (str): The directory containing the data files.
        part (str): The part of the dataset (e.g., 'train', 'val', 'test').

    Returns:
        dict: A dictionary mapping SMILES strings to their corresponding
              atom map numbers and positions.
    """
    conf_dir = os.path.join(data_dir, f'conf_{part}')
    smiles2conf = {}
    with open(os.path.join(data_dir, f'amlist_{part}.txt')) as Fin:
        all_info = Fin.readlines()
        for lin in tqdm(all_info):
            if len(lin) < 4:
                continue
            smiles, charge, name, amlist = lin.strip().split('\t')
            if os.path.exists(os.path.join(conf_dir, name)):
                amlist = eval(amlist)
                mol3d = aseread(os.path.join(conf_dir, name))
                positions = mol3d.get_positions()
                am2pos = {v: positions[k].tolist() for k, v in amlist}
                # use to list to save the space for store
                smiles2conf[smiles] = am2pos



def load_mol_with_conf(smiles:str, smiles2conf:dict, conf_dir:str, addHs=False):
    '''
    将smiles转化为mol对象，并且加载对应的conf：返回mol对象、构象、smiles part mask
    Args:
        smiles: smiles
        smiles2conf: smiles2conf dict
        conf_dir: conf_dir
    Returns:
        mol: mol对象
        conf: 构象
        smiles_part_mask: N * N mask（表示smiles中不同的分子）
    '''

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    if addHs:
        mol =Chem.AddHs(mol)
        extra_Hs = []

    all_atom_numbers = 0
    single_part_atom_numbers = []
    all_confs = []
    # debug:
    # org_atom_list = np.array([x.GetAtomicNum() for x in mol.GetAtoms()])
    all_atoms = []
    single_parts = smiles.split('.')
    for smiles in single_parts:
        cano_smiles = clear_map_number(smiles)
        if cano_smiles not in smiles2conf:
            assert False, f'no conf for {smiles}'
        cano_mol_noHs = Chem.MolFromSmiles(cano_smiles)
        conf_path = path_join(conf_dir, smiles2conf[cano_smiles])
        conf = aseread(conf_path)
        # all_atoms.append(conf.get_atomic_numbers())
        conf = conf.get_positions()
        atom_numbers = len(cano_mol_noHs.GetAtoms())

        if not addHs:
            all_confs.append(conf[:atom_numbers])
            all_atom_numbers += atom_numbers
            single_part_atom_numbers.append(atom_numbers)
            # all_atoms[-1] = all_atoms[-1][:atom_numbers]
        else:
            all_confs.append(conf)
            all_atom_numbers += len(conf)
            single_part_atom_numbers.append(len(conf))
            extra_Hs.append(len(conf) - atom_numbers)

    # all_atoms = np.concatenate(all_atoms)

    mol_part_mask = np.zeros((all_atom_numbers, all_atom_numbers), dtype=np.float32)
    for i in range(len(single_part_atom_numbers)):
        start = sum(single_part_atom_numbers[:i])
        end = start + single_part_atom_numbers[i]
        mol_part_mask[start:end, start:end] = 1.0

    all_confs = np.concatenate(all_confs, axis=0)
    if addHs:
        # rerank conf (mol1, mol1_Hs, ..., mol_n, mol_n_Hs -> mol1, mol2, ..., mol1_Hs, ..., mol_n_Hs) and mask
        # that is: move the extra Hs to the end of the conf
        noHs_idx = []
        cur_idx = 0
        for total_num, Hs_num in zip(single_part_atom_numbers, extra_Hs):
            noHs_idx.append(cur_idx + np.arange(total_num-Hs_num))
            cur_idx += total_num
        noHs_idx = np.concatenate(noHs_idx)
        Hs_idx = np.setdiff1d(np.arange(all_atom_numbers), noHs_idx)
        rerank_idx = np.concatenate([noHs_idx, Hs_idx])
        # rerank conf
        all_confs = all_confs[rerank_idx]
        # rerank mask
        mol_part_mask = mol_part_mask[rerank_idx][:, rerank_idx]
        # debug
        # all_atoms = all_atoms[rerank_idx]

    # assert np.array_equal(all_atoms, org_atom_list), \
    #     f'atomic number not align: {smiles} {conf_dir} {org_atom_list} {all_atoms}'

    return mol, all_confs, mol_part_mask


if __name__ == '__main__':
    from utils import clear_map_number
    import json, pandas as pd
    parts = ['test', 'train', 'val']

    for part in parts:
        smi2conf_path = f'data/smi2conf_{part}.json'
        conf_dir = f'data/conf_{part}'
        with open(smi2conf_path, 'r') as f:
            smi2conf = json.load(f)

        df = pd.read_csv(f'data/{part}.csv')
        rxns = df['mapped_rxn'].to_list()

        for rxn in rxns:
            # rxn = ['[O:1]=[C:2](/[N:3]=[CH:4]/[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1)[c:18]1[cH:19][cH:20][cH:21][cH:22][cH:23]1.[SH:5][c:6]1[cH:7][cH:8][cH:9][cH:10][cH:11]1>>[O:1]=[C:2]([NH:3][CH:4]([S:5][c:6]1[cH:7][cH:8][cH:9][cH:10][cH:11]1)[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1)[c:18]1[cH:19][cH:20][cH:21][cH:22][cH:23]1']
            r, p = rxn.split('>>')

            for smi in [r, p]:
                mol, conf, mask = load_mol_with_conf(smi, smi2conf, conf_dir, addHs=True)
                # print(len(mask), len(conf), len(mol.GetAtoms()))
                assert len(mask) == len(conf) == len(mol.GetAtoms()), f'len not equal: {len(mask)} {len(conf)} {len(mol.GetAtoms())}'
                print(conf, mask)
                exit(0)