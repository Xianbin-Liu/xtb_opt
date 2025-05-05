import json
from tqdm import tqdm
import os
import rdkit
from rdkit import Chem
import argparse
import multiprocessing
from rdkit.Chem import AllChem
import pandas
from utils import cano, clear_map_number, cano_readd


# Rkey = 'reactants>reagents>production'
Rkey = 'mapped_rxn'

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

def gene_conf(args):
    sms, tlock, globalid, out_dir, split = args
    mol = Chem.MolFromSmiles(sms)
    if mol is None:
        return False, {'name': sms}
    mhs = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20000819
    params.useRandomCoords = True
    params.maxAttempts = 500
    succ = AllChem.EmbedMolecule(mhs, params)

    if succ == -1:
        return False, {'name': sms}

    with tlock:
        globalid.value += 1
        cval = globalid.value

    if split > 0:
        tdx = (cval - 1) // split
        mid_path = f'{tdx * split + 1}-{tdx * (split + 1)}'
        out_dir = os.path.join(out_dir, mid_path)

    Chem.MolToXYZFile(mhs, os.path.join(out_dir, f'conf_{cval}.xyz'))
    # am_list = [
    #     (idx, x.GetAtomMapNum()) for idx, x in enumerate(mhs.GetAtoms())
    #     if x.HasProp('molAtomMapNumber')
    # ]

    total_charge = sum(x.GetFormalCharge() for x in mol.GetAtoms())

    return True, {
        'charge': total_charge,
        'conf_name': f'conf_{cval}.xyz', 'name': sms
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Parser')
    parser.add_argument(
        '--data_dir', type=str, default='data', #required=True,
        help='the dir containing data path'
    )

    parser.add_argument(
        '--output_dir', type=str, default='data', #required=True,
        help='the dir to output data'
    )
    parser.add_argument(
        '--nproc', type=int, default=1,
        help='the number of process'
    )
    parser.add_argument(
        '--split', type=int, default=-1,
        help='the split folder to store conf, ' +
        'non-positive for not splitting'
    )
    args = parser.parse_args()

    Man = multiprocessing.Manager()
    tlock = Man.Lock()
    Vl = Man.Value('i', 0)

    for part in ['val', 'train', 'test']:
        yargs, xargs, outlines, amaplist, failed = set(), [], [], [], []
        Pl, added = multiprocessing.Pool(args.nproc), set()
        data_path = os.path.join(args.data_dir, f'{part}.csv')
        df = pandas.read_csv(data_path)
        co_path = os.path.join(args.output_dir, f'conf_{part}')
        os.makedirs(co_path, exist_ok=True)
        all_rxns = df[Rkey].tolist()

        # extract all mols
        no_am_cano_smis = set()
        for rxn in all_rxns:
            reacs, prods = rxn.split('>>')
            reacs_list = clear_map_number(reacs).split('.')
            prods_list = clear_map_number(prods).split('.')
            for reac in reacs_list + prods_list:
                no_am_cano_smis.add(cano(reac))

        yargs = [(x, tlock, Vl, co_path, args.split) for x in (no_am_cano_smis)]
        supp_iter = Pl.imap_unordered(gene_conf, yargs)


        smi2conf = {}

        for succ, info in tqdm(supp_iter, total=len(yargs)):
            if not succ:
                continue
            amaplist.append((
                info['name'], str(info['charge']),
                info['conf_name']
            ))
            smi2conf[info['name']] = info['conf_name']

        with tlock:
            Vl.value = 0

        Pl.close()
        Pl.join()

        p_path = os.path.join(args.output_dir, f'processed_{part}.csv')
        am_path = os.path.join(args.output_dir, f'amlist_{part}.txt')
        f_path = os.path.join(args.output_dir, f'failed_{part}.csv')
        smi2conf_path = os.path.join(args.output_dir, f'smi2conf_{part}.json')
        with open(smi2conf_path, 'w') as F:
            json.dump(smi2conf, F, indent=4)

        with open(am_path, 'w') as F:
            for line in amaplist:
                F.write('{}\n'.format('\t'.join(line)))

        if len(outlines) > 0:
            df_out = pandas.DataFrame(outlines)
            df_out.to_csv(p_path, index=False)

        if len(failed) > 0:
            faild_out = pandas.DataFrame(failed)
            faild_out.to_csv(f_path, index=False)
