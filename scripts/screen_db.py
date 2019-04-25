#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 12.05.16
# license         : BSD-3
# ==============================================================================

import os
import sys
import time
import argparse
import marshal
import sqlite3 as lite
from collections import defaultdict
from pmapper.pharmacophore import Pharmacophore


def compare_fp(query_fp, fp):
    return (query_fp & fp) == query_fp


def main(db_fname, query_fname, out_fname):
    if out_fname is not None:
        out_f = open(out_fname, 'wt')

    # load query
    q = Pharmacophore(cached=True)
    q.load_from_pma(query_fname)
    q_fp = q.get_fp()

    conn = lite.connect(db_fname)
    cur = conn.cursor()

    cur.execute("SELECT bin_step FROM settings")
    db_bin_step = cur.fetchone()[0]

    # check bin steps
    if q.get_bin_step() != db_bin_step:
        sys.stderr.write('Model has a different bin step from compounds in database. '
                         'It would be skipped.\n')
        raise Exception('Incompatible bin step')

    cur.execute("SELECT DISTINCT(mol_name) FROM conformers")
    mol_names = [i[0] for i in cur.fetchall()]

    for mol_name in mol_names:

        # load fp for all conformers
        cur.execute("SELECT conf_id, fp FROM conformers WHERE mol_name = ?", (mol_name,))
        data = cur.fetchall()

        conf_ids = []
        for conf_id, fp in data:
            if compare_fp(q_fp, marshal.loads(fp)):
                conf_ids.append(conf_id)

        if conf_ids:
            # load pharmacophores for selected conformers
            sql = "SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN (%s)" % \
                  ','.join(['?'] * len(conf_ids))
            cur.execute(sql, conf_ids)
            res = cur.fetchall()
            confs = defaultdict(list)
            for r in res:
                confs[r[0]].append((r[1], tuple(r[2:])))
            for conf_id, feature_coords in confs.items():

                p = Pharmacophore(bin_step=db_bin_step, cached=True)
                p.load_from_feature_coords(feature_coords)
                res = p.fit_model(q)
                if res:
                    res_string = '\t'.join(map(str, (mol_name, conf_id, ','.join(map(str, res))))) + '\n'
                    if out_fname is not None:
                        out_f.write(res_string)
                    else:
                        sys.stdout.write(res_string)
                    break  # match the first conformer


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries. '
                                                 'Output is printed out in STDOUT.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database', metavar='input.db', required=True,
                        help='SQLite input file with compounds to screen.')
    parser.add_argument('-q', '--query', metavar='model.pma', required=True,
                        help='pharmacophore models. Several models can be specified. Model format is recognized by '
                             'file extension. Currently pma (pmapper) and pml (LigandScout) formats are supported.')
    parser.add_argument('-o', '--output', default=None,
                        help='path to output text file with screening results.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "database": db_fname = v
        if o == "query": query_fname = v
        if o == "output": out_fname = v

    main(db_fname=db_fname,
         query_fname=query_fname,
         out_fname=out_fname)