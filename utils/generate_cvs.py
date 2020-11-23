#!/usr/bin/env python
import argparse
from parmed.amber import AmberParm, AmberMask


def get_atom_index(top, mask):
    return next(AmberMask(top, mask).Selected()) + 1


def write_cv(top, raw_cv_def):
    cv_def = raw_cv_def.split()
    cv_type = cv_def[0]
    indices = [str(get_atom_index(top, mask)) for mask in cv_def[1:]]
    print('\n$COLVAR')
    print('COLVAR_type = "{}"'.format(cv_type))
    print('atoms = {}'.format(', '.join(indices)))
    print('$END')


def main():
    parser = argparse.ArgumentParser(description="Create CVs file from amber masks")
    parser.add_argument("top", help="Topology file")
    parser.add_argument("defs", help="CV definitions")
    args = parser.parse_args()

    top = AmberParm(args.top)
    with open(args.defs, 'r') as f:
        cv_defs = f.readlines()

    print('{}'.format(len(cv_defs)))
    for cv_def in cv_defs:
        write_cv(top, cv_def)


if __name__ == '__main__':
    main()
