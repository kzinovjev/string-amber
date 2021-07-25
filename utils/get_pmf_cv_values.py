#!/usr/bin/env python

import os
import re


REORDER_TRJ_OUTPUT_DIR = 'reorder_trj_results'
IN_FILE = '1.in'
REORDER_TRJ_FORMAT = 'nc'


def get_ntwx():
    with open(IN_FILE) as f:
        ntwx_regex = re.compile(r".*(^|[, \t]+)ntwx[ \t]*=[ \t]*([0-9]*)[\D]*")
        for line in f:
            match = ntwx_regex.match(line.strip())
            if match:
                return int(match.groups()[1])


def get_results_dir():
    with open('STRING') as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith('dir'):
                return line.split()[-1].strip('"')
    return '.'


def get_nnodes():
    os.listdir(os.path.join(REORDER_TRJ_OUTPUT_DIR, 'pmf'))
    paths = filter(lambda name: name.endswith(REORDER_TRJ_FORMAT),
                   os.listdir(os.path.join(REORDER_TRJ_OUTPUT_DIR, 'pmf')))
    return len(list(paths))


def read_pmf_dat(filename, ntwx):
    with open(filename, 'r') as f:
        return f.readlines()[ntwx+2::ntwx]


def get_out_path(node):
    return os.path.join(REORDER_TRJ_OUTPUT_DIR,
                        'pmf',
                        '{}.dat'.format(node))


def main():
    for node in range(1, get_nnodes() + 1):
        print("Extracting from {}_final.dat".format(node))
        dat_filename = os.path.join(get_results_dir(), '{}_final.dat'.format(node))

        with open(get_out_path(node), 'w') as f:
            for line in read_pmf_dat(dat_filename, get_ntwx()):
                f.write(line)


if __name__ == '__main__':
    main()
