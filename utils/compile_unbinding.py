#!/usr/bin/env python
import re
import argparse
from parmed.amber import AmberParm, AmberMask


CVS = (('r_r',         'MBOND',     ('P1', 'L1')),
       ('r_theta',     'MANGLE',    ('P1', 'L1', 'L2')),
       ('r_phi',       'MDIHEDRAL', ('P1', 'L1', 'L2', 'L3')),
       ('omega_theta', 'MANGLE',    ('P2', 'P1', 'L1')),
       ('omega_phi',   'MDIHEDRAL', ('P2', 'P1', 'L1', 'L2')),
       ('omega_psi',   'MDIHEDRAL', ('P3', 'P2', 'P1', 'L1')))


def read_mask(line):
    center_re = re.compile(r"([^:]*):(.*)")
    return [entry.strip() for entry in center_re.match(line).groups()]


def read_masks(filename):
    with open(filename, 'r') as f:
        return {name: mask for name, mask in map(read_mask, f)}


def get_center(top, mask):
    return [str(i + 1) for i in AmberMask(top, mask).Selected()]


def write_centers(filename, centers):
    with open(filename, 'w') as f:
        for center in centers:
            f.write(str(len(center)) + '\n')
            f.write(' '.join(center) + '\n\n')


def filter_centers(centers, center_list):
    return [centers[center] for center in center_list]


def compile_cv(cvs_file, centers, name, kind, center_list):
    cv_filename = '{}.cv'.format(name)
    cvs_file.write('$COLVAR\n')
    cvs_file.write('COLVAR_type = "{}"\n'.format(kind))
    cvs_file.write('atoms_file = "{}"\n'.format(cv_filename))
    cvs_file.write('$END\n\n')
    write_centers(cv_filename, filter_centers(centers, center_list))


def main():
    parser = argparse.ArgumentParser(description="Create CVs file for unbinding")
    parser.add_argument("top", help="Topology file")
    parser.add_argument("masks", help="Masks file (with masks for 6 centers (P1..L3")
    args = parser.parse_args()

    top = AmberParm(args.top)
    masks = read_masks(args.masks)
    centers = {name: get_center(top, mask) for name, mask in masks.items()}

    with open('CVs', 'w') as f:
        f.write('6\n\n')
        for cv in CVS:
            compile_cv(f, centers, *cv)


if __name__ == '__main__':
    main()
