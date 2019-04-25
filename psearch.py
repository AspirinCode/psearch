#!/usr/bin/env python3

import os
import sys
import argparse

from scripts import gen_pharm_models
from scripts import select_training_set_rdkit
from scripts import screen_db
import external_statistics


def _make_dir(new_path):
    paths = [os.path.split(os.path.split(new_path)[0])[0], os.path.split(new_path)[0], new_path]
    for path in paths:
        if not os.path.exists(path):
            os.mkdir(path)


def calc(mol_act, mol_inact, in_adb, in_indb, files_at, files_int, tol, lower):
    path_pma, ids = gen_pharm_models.main(in_adb=in_adb,
                                          in_indb=in_indb,
                                          act_trainset=files_at,
                                          inact_trainset=files_int,
                                          out_pma=None,
                                          tolerance=tol,
                                          lower=lower)

    if not path_pma:
        return 'Error: {}'.format(os.path.dirname(files_at))
    elif ids == 0:
        return 'Error: {}'.format(path_pma)

    path_screen = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_adb)))[0], 'screen',
                               '{}_ph{}'.format(os.path.basename(files_at).split('.')[0].split('_')[1], ids))
    _make_dir(path_screen)

    for query_fname in os.listdir(path_pma):
        out_f_screen_act = os.path.join(path_screen, 'screen_{}_{}.txt'.format(
            os.path.basename(in_adb).split('.')[0], query_fname.split('.')[0]))
        out_f_screen_inact = os.path.join(path_screen, 'screen_{}_{}.txt'.format(
            os.path.basename(in_indb).split('.')[0], query_fname.split('.')[0]))

        screen_db.main(db_fname=in_adb, query_fname=os.path.join(path_pma, query_fname), out_fname=out_f_screen_act)
        screen_db.main(db_fname=in_indb, query_fname=os.path.join(path_pma, query_fname), out_fname=out_f_screen_inact)

    out_external_stat = os.path.join(os.path.split(os.path.dirname(os.path.abspath(in_adb)))[0],
                                     'external_statistics{}.txt'.format(os.path.split(path_pma)[1]))

    external_statistics.main(mol_act=mol_act,
                             mol_inact=mol_inact,
                             ts_act=files_at,
                             ts_inact=files_int,
                             path_to_pma=path_pma,
                             path_to_screen=path_screen,
                             out_external=out_external_stat)


def main(mol_act, mol_inact, in_adb, in_indb, mode_train_set,
         tol, lower, fdef_fname, treshold_clust, clust_size, max_nact_trainset):
    if 1 in mode_train_set:
        list_ts_1 = select_training_set_rdkit.main(in_fname_act=mol_act,
                                                   in_fname_inact=mol_inact,
                                                   fdef_fname=fdef_fname,
                                                   make_clust=False,
                                                   fcfp4=fcfp4,
                                                   clust_stat=None,
                                                   treshold_clust=treshold_clust,
                                                   clust_size=clust_size,
                                                   max_nact_trainset=max_nact_trainset)
    else:
        list_ts_1 = []

    if 2 in mode_train_set:
        list_ts_2 = select_training_set_rdkit.main(in_fname_act=mol_act,
                                                   in_fname_inact=mol_inact,
                                                   fdef_fname=fdef_fname,
                                                   make_clust=True,
                                                   fcfp4=fcfp4,
                                                   clust_stat=open(os.path.join(
                                                       os.path.dirname(os.path.abspath(mol_act)),
                                                       'cluster_stat_thr{}.txt'.format(treshold_clust)), 'w'),
                                                   treshold_clust=treshold_clust,
                                                   clust_size=clust_size,
                                                   max_nact_trainset=max_nact_trainset)
    else:
        list_ts_2 = []

    list_ts = list_ts_1 + list_ts_2

    for files_at, files_int in list_ts:
        res = calc(mol_act, mol_inact,
                   in_adb, in_indb,
                   files_at, files_int,
                   tol, lower)
        if res:
            sys.stderr.write(res)
        continue


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-mol_act', '--active_mol', metavar='active.smi', required=True,
                        help='.smi file with active molecules.')
    parser.add_argument('-mol_inact', '--inactive_mol', metavar='inactive.smi', required=True,
                        help='.smi file with inactive molecules.')
    parser.add_argument('-iadb', '--input_active_db', metavar='input.db', required=True,
                        help='SQLite DB with active pharmacophores (feature coordinates).')
    parser.add_argument('-iindb', '--input_inactive_db', metavar='input.db', required=True,
                        help='SQLite DB with inactive pharmacophores (feature coordinates).')
    parser.add_argument('-ts', '--mode_train_set', nargs='+', default=[1, 2], type=int,
                        help='1 - form a training set by Stategy 1,'
                             '2 - form a training set by Stategy 2,'
                             '1 2 - form a training sets by Stategy 1 and Stategy 2,')
    parser.add_argument('-tol', '--tolerance', metavar='VALUE', default=0,
                        help='tolerance used for calculation of a stereoconfiguration sign.')
    parser.add_argument('-l', '--lower', metavar='3', type=int, default=3,
                        help='lower number of features used for generation of subpharmacophores.')
    parser.add_argument('--fcfp4', action='store_true', default=False,
                        help='if set FCFP4 fingerprints will be used for compound selection, '
                             'otherwise pharmacophore fingerprints will be used based on feature '
                             'definitions provided by --rdkit_fdef argument.')
    parser.add_argument('-t', '--treshold_clust', default=0.4,
                        help='treshold for сlustering data by Butina algorithm')
    parser.add_argument('-clz', '--clust_size', default=5,
                        help='minimum cluster size from extract centroinds for training set')
    parser.add_argument('-ma', '--max_act_ts', default=5,
                        help='maximum number of active compounds for training set')
    parser.add_argument('--fdef', metavar='smarts.fdef', required=True,
                        help='fdef-file with pharmacophore feature definition.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "active_mol": active_mol = v
        if o == "inactive_mol": inactive_mol = v
        if o == "input_active_db": in_adb = v
        if o == "input_inactive_db": in_indb = v
        if o == "mode_train_set": mode_train_set = v
        if o == "tolerance": tol = int(v)
        if o == "lower": lower = int(v)
        if o == "input_active_mol": mol_act = v
        if o == "input_inactive_mol": mol_inact = v
        if o == "fcfp4": fcfp4 = v
        if o == "treshold_clust": treshold_clust = float(v)
        if o == "clust_size": clust_size = int(v)
        if o == "max_act_ts": max_nact_trainset = int(v)
        if o == "fdef": fdef_fname = v

    main(mol_act=active_mol,
         mol_inact=inactive_mol,
         in_adb=in_adb,
         in_indb=in_indb,
         mode_train_set=mode_train_set,
         tol=tol,
         lower=lower,
         treshold_clust=treshold_clust,
         clust_size=clust_size,
         max_nact_trainset=max_nact_trainset,
         fdef_fname=fdef_fname)