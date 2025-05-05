import json
from tqdm import tqdm
import os
import rdkit
from rdkit import Chem
import argparse
import multiprocessing
from rdkit.Chem import AllChem
import pandas

# Rkey = 'reactants>reagents>production'
Rkey = 'mapped_rxn'


def cano(x):
    return Chem.MolToSmiles(Chem.MolFromSmiles(x))


def clear_map_number(x):
    mol = Chem.MolFromSmiles(x)
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atom.ClearProp('molAtomMapNumber')
    tmp_x = Chem.MolToSmiles(mol)
    return Chem.MolToSmiles(Chem.MolFromSmiles(tmp_x))


def cano_readd(x):
    x = clear_map_number(x)
    mol = Chem.MolFromSmiles(x)
    for idx, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(idx + 1)
    return Chem.MolToSmiles(mol)


def remove_same_useless(rxn):
    reac, prod = rxn.split('>>')
    reac = [cano(x) for x in reac.split('.')]
    prod = [cano(x) for x in prod.split('.')]

    cnt_x, cnt_y = {}, {}
    for x in reac:
        cnt_x[x] = cnt_x.get(x, 0) + 1

    for x in prod:
        cnt_y[x] = cnt_y.get(x, 0) + 1

    real_prod = []
    for x in prod:
        if cnt_x.get(x, 0) == 0:
            real_prod.append(x)
        else:
            cnt_x[x] -= 1

    real_reac = [k for k, v in cnt_x.items() if v > 0]

    am_prod = [
        x.GetAtomMapNum() for y in real_prod
        for x in Chem.MolFromSmiles(y).GetAtoms()
    ]

    prelen_amprod = len(am_prod)
    am_prod = set(am_prod)

    if prelen_amprod != len(am_prod):
        return False, "Duplicated Atom Mapping in Prod"
    if 0 in am_prod:
        return False, 'Product Atom not mapped'

    am_reac = set()

    for x in real_reac:
        this_mol = Chem.MolFromSmiles(x)
        this_am = [
            t.GetAtomMapNum() for t in this_mol.GetAtoms()
            if t.HasProp('molAtomMapNumber')
        ]

        prelen = len(this_am)
        this_am = set(this_am)
        if prelen != len(this_am) or len(this_am & am_reac) != 0:
            return False, 'Duplicated Atom Mapping in Reac'
        am_reac.update(this_am)

    if len(am_prod - am_reac) != 0:
        return False, 'Prod Atom Mapping not in Reac Found'

    if len(real_prod) != 1:
        return False, 'Multiple Product Given'

    max_am_reac = max(am_reac)
    amcmp_real_reac = []
    for x in real_reac:
        this_mol = Chem.MolFromSmiles(x)
        for t in this_mol.GetAtoms():
            if not t.HasProp('molAtomMapNumber'):
                t.SetAtomMapNum(max_am_reac + 1)
                max_am_reac += 1
        amcmp_real_reac.append(Chem.MolToSmiles(this_mol))

    prod_mol = Chem.MolFromSmiles(real_prod[0])
    cano_prod = Chem.MolToSmiles(prod_mol, ignoreAtomMapNumbers=True)

    def remap_am(x, remap):
        mol = Chem.MolFromSmiles(x)
        for x in mol.GetAtoms():
            if x.HasProp('molAtomMapNumber'):
                tmap = x.GetAtomMapNum()
                remap[tmap] = remap.get(tmap, len(remap) + 1)
                x.SetAtomMapNum(remap[tmap])
        return Chem.MolToSmiles(mol)

    ## TODO：问题：如果抵消的”same part“的am是否还正确：【0，1，2，3，4】 -》 【0，2】；通过remap-am就会重新在2开始填补
    xmol, remap = Chem.MolFromSmiles(cano_prod), {}
    for atom in xmol.GetAtoms():
        tmap = atom.GetAtomMapNum()
        remap[tmap] = len(remap) + 1
        atom.SetAtomMapNum(remap[tmap])
    prod_out = Chem.MolToSmiles(xmol)
    amcmp_real_reac = [remap_am(x, remap) for x in amcmp_real_reac]

    return True, ('.'.join(amcmp_real_reac), prod_out)


def process_one(args):
    x, classid, idx, tlock, globalid, out_dir, split = args
    reac, prod = x.split('>>')
    if Chem.MolFromSmiles(reac) is None or \
            Chem.MolFromSmiles(prod) is None:
        return False, {
            Rkey: x, 'class': classid, 'id': idx,
            'msg': "Invalid SMILES"
        }

    Succ, Msg = remove_same_useless(x)
    if not Succ:
        return False, {Rkey: x, 'class': classid, 'id': idx, 'msg': Msg}

    reac, prod = Msg

    Pmol = Chem.MolFromSmiles(prod)
    if Pmol.GetNumAtoms() == 1:
        return False, {
            Rkey: x, 'class': classid, 'id': idx,
            'msg': 'Single Atom Product'
        }

    mhs = Chem.AddHs(Pmol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20000819
    params.useRandomCoords = True
    params.maxAttempts = 500
    succ = AllChem.EmbedMolecule(mhs, params)

    if succ == -1:
        return False, {
            Rkey: x, 'class': classid, 'id': idx,
            'msg': "No Conformation"
        }
    with tlock:
        globalid.value += 1
        cval = globalid.value

    if split > 0:
        tdx = (cval - 1) // split
        mid_path = f'{tdx * split + 1}-{tdx * (split + 1)}'
        out_dir = os.path.join(out_dir, mid_path)

    Chem.MolToXYZFile(mhs, os.path.join(out_dir, f'conf_{cval}.xyz'))
    am_list = [
        (idx, x.GetAtomMapNum()) for idx, x in enumerate(mhs.GetAtoms())
        if x.HasProp('molAtomMapNumber')
    ]

    total_charge = sum(x.GetFormalCharge() for x in Pmol.GetAtoms())

    return True, {
        Rkey: f'{reac}>>{prod}', 'class': classid, 'id': idx,
        'am_list': am_list, 'charge': total_charge,
        'conf_name': f'conf_{cval}.xyz', 'product': prod
    }


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
    am_list = [
        (idx, x.GetAtomMapNum()) for idx, x in enumerate(mhs.GetAtoms())
        if x.HasProp('molAtomMapNumber')
    ]

    total_charge = sum(x.GetFormalCharge() for x in mol.GetAtoms())

    return True, {
        'am_list': am_list, 'charge': total_charge,
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
        for idx, row in df.iterrows():
            xargs.append((
                row[Rkey], row.get('class', -1), row['id'],
                tlock, Vl, co_path, args.split
            ))

        res_iter = Pl.imap_unordered(process_one, xargs)
        for succ, info in tqdm(res_iter, total=len(xargs)):
            if not succ:
                failed.append(info)
            else:
                outlines.append({k: info[k] for k in [Rkey, 'class', 'id']})
                amaplist.append((
                    info['product'], str(info['charge']),
                    info['conf_name'], json.dumps(info['am_list'])
                ))
                reac_list = info[Rkey].split('>>')[0].split('.')
                yargs.update(cano_readd(x) for x in reac_list)
                added.add(info[Rkey].split('>>')[1])

        yargs = [(x, tlock, Vl, co_path, args.split) for x in (yargs - added)]
        supp_iter = Pl.imap_unordered(gene_conf, yargs)
        for succ, info in tqdm(supp_iter, total=len(yargs)):
            if not succ:
                continue
            amaplist.append((
                info['name'], str(info['charge']),
                info['conf_name'], json.dumps(info['am_list'])
            ))

        with tlock:
            Vl.value = 0

        Pl.close()
        Pl.join()

        p_path = os.path.join(args.output_dir, f'processed_{part}.csv')
        am_path = os.path.join(args.output_dir, f'amlist_{part}.txt')
        f_path = os.path.join(args.output_dir, f'failed_{part}.csv')

        if len(outlines) > 0:
            df_out = pandas.DataFrame(outlines)
            df_out.to_csv(p_path, index=False)
            with open(am_path, 'w') as F:
                for line in amaplist:
                    F.write('{}\n'.format('\t'.join(line)))

        if len(failed) > 0:
            faild_out = pandas.DataFrame(failed)
            faild_out.to_csv(f_path, index=False)
