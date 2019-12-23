#!/usr/bin/env python

import mdtraj as md
import argparse
import os
from functools import reduce


def sum_rsts(rsts):
    return reduce(lambda a, b: a + b, rsts)


def load_rst(filename, top):
    print("Loading {}".format(filename))
    try:
        return md.load_ncrestrt(filename, top=top)
    except (IOError, TypeError):
        pass
    try:
        return md.load_restrt(filename, top=top)
    except (IOError, TypeError):
        print('Trajectory format not recognized. Exiting.')
        exit()


def get_rst_indexes(rex_file):
    with open(rex_file, 'r') as f:
        rex_data = f.readlines()[-1].split()[1:]
        return [rex_data.index(str(i + 1)) + 1 for i in range(len(rex_data))]


def get_rst(index, ext, top):
    return load_rst('{}.{}'.format(index, ext), top)


def get_rsts(indexes, ext, top):
    return [get_rst(i, ext, top) for i in indexes]


def main():

    parser = argparse.ArgumentParser(
        description=""
        "Creates a trajectory of the last frames of all (numbered)"
        "trajectories in the current directory.\n",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("prmtop", help="Amber topology file")
    parser.add_argument("rex", help="REX plot file (plot.REX or plot_final.REX)")
    parser.add_argument("output", help="Name of the output trajectory")
    parser.add_argument("--ext", help="extension of the rst files", 
                        default='rst')
    args = parser.parse_args()

    rsts = get_rsts(get_rst_indexes(args.rex), args.ext, args.prmtop)
    sum_rsts(rsts).save(args.output)


if __name__ == '__main__':
    main()
