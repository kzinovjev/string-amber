#!/usr/bin/env python

import argparse
import molmod as mm
import numpy as np
import mdtraj as md


def eval_pplane(atoms_xyz):
    return mm.opbend_dist(atoms_xyz)[0] / mm.angstrom


CV_KIND_DICT = {
    'BOND': lambda atoms_xyz: mm.bond_length(atoms_xyz)[0] / mm.angstrom,
    'ANGLE': lambda atoms_xyz: mm.bend_angle(atoms_xyz)[0] / mm.deg,
    'DIHEDRAL': lambda atoms_xyz: mm.dihed_angle(atoms_xyz)[0] / mm.deg,
    'PPLANE': lambda atoms_xyz: eval_pplane(np.roll(atoms_xyz, -3))
}


def eval_cv(xyz, atoms, callback):
    return callback([xyz[atom-1]*mm.nanometer for atom in atoms])


def get_cv_value(xyz, atoms, kind):
    return eval_cv(xyz, atoms, CV_KIND_DICT[kind])


def extract_value(lines, key):
    """
    :param lines: list of lines
    :param key: substring that the line with the value must contain
    :return: everything after '=' from a first line in lines that contains key
    """
    return filter(lambda line: key in line, lines)[0].split('= ')[1].strip('\n')


def parse_cv(raw_cv):
    """
    :param raw_cv: raw text lines defining a single CV
    :return: a callable that accepts xyz and returns the CV value
    """
    kind = extract_value(raw_cv, 'COLVAR_type').strip('"')
    atoms = [int(idx) for idx in extract_value(raw_cv, 'atoms').split(',')]
    return lambda xyz: eval_cv(xyz, atoms, CV_KIND_DICT[kind])


def lines_from_file(filename):
    with open(filename, 'r') as f:
        return f.readlines()


def read_raw_cvs(filename):
    """
    :param filename: CVs filename
    :return: list of subsets of lines from filename corresponding to each CV
    """
    lines = lines_from_file(filename)
    start_lines = [idx for idx, line in enumerate(lines) if '$COLVAR' in line]
    ranges = zip(start_lines, start_lines[1:] + [len(lines),])
    return [lines[start:end] for start, end in ranges]


def read_cvs(filename):
    return map(parse_cv, read_raw_cvs(filename))


def load_trj(filename, top):
    try:
        return md.load(filename, top=top)
    except IOError:
        pass
    try:
        return md.load_netcdf(filename, top=top)
    except IOError:
        pass
    try:
        return md.load_mdcrd(filename, top=top)
    except IOError:
        print('Trajectory format not recognized. Exiting.')
        exit()


def xyz_to_cv_values(xyz, cvs):
    return tuple(cv(xyz) for cv in cvs)


def print_cv_values(cv_values):
    print('{:13.5e}' * len(cv_values)).format(*cv_values)


def main():
    parser = argparse.ArgumentParser(description="Extract CV values from a trajectory")
    parser.add_argument("cvs_file", help="CVs file")
    parser.add_argument("trj_file", help="Trajectory file (any format supported by mdtraj)")
    parser.add_argument("top", help="Topology file")
    parser.add_argument("first", type=int, help="Start frame (default: 1)", default=1, nargs='?')
    parser.add_argument("last", type=int, help="End frame (default: <last frame>)", nargs='?')
    args = parser.parse_args()

    cvs = read_cvs(args.cvs_file)
    trj = load_trj(args.trj_file, args.top)[args.first-1:args.last]

    cv_values = np.array([xyz_to_cv_values(frame, cvs) for frame in trj.xyz])

    for frame_cvs in cv_values:
        print_cv_values(frame_cvs)
    print('\nAVERAGE:')
    print_cv_values(np.mean(cv_values, axis=0))


if __name__ == '__main__':
    main()
