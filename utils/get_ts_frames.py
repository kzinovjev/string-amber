#!/usr/bin/env python

import os
import argparse
import numpy as np
import re


REORDER_TRJ_OUTPUT_DIR = 'reorder_trj_results'
IN_FILE = '1.in'
CPPTRAJ_IN_NAME = '.reorder_cpptraj_tmp.in'


def assert_file_exist(filename):
    if not os.path.exists(filename):
        print('{} not found. Exiting.'.format(filename))
        exit()


def get_ntwx():
    assert_file_exist(IN_FILE)
    with open(IN_FILE) as f:
        ntwx_regex = re.compile(r".*=(\D*)([0-9]*)")
        for line in f:
            if line.strip().startswith('ntwx'):
                return int(ntwx_regex.match(line.strip()).groups()[1])


def get_results_dir():
    assert_file_exist('STRING')
    with open('STRING') as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith('dir'):
                return line.split()[-1].strip('"')
    return '.'


def get_nnodes():
    paths = filter(lambda name: name.endswith('nc'),
                   os.listdir(os.path.join(REORDER_TRJ_OUTPUT_DIR, 'pmf')))
    return len(list(paths))


def interpolate_max(x, y):
    x_matrix = np.array([x**2, x, np.ones(3)]).T
    c = np.linalg.inv(x_matrix).dot(y)
    return -c[1]/(c[0]*2)


def get_ts(pmf_filename):
    pmf = np.loadtxt(pmf_filename, skiprows=2, usecols=(0, 5)).T
    max_i = np.argmax(pmf, 1)[1]
    x = pmf[0][max_i-1:max_i+2]
    y = pmf[1][max_i-1:max_i+2]
    return interpolate_max(x, y)


def read_pmf_dat(filename, ntwx):
    return np.loadtxt(filename, skiprows=3, usecols=0)[ntwx-1::ntwx]


def get_ts_frames(s_list, ts, thr):
    return [i + 1 for i, s in enumerate(s_list) if np.abs(s-ts) < thr]


def parse_pmf_dat(results_dir, ntwx, ts, thr):
    ts_frames = {}
    for node in range(1, get_nnodes() + 1):
        dat_filename = os.path.join(results_dir, '{}_final.dat'.format(node))
        s = read_pmf_dat(dat_filename, ntwx)
        ts_frames[node] = get_ts_frames(s, ts, thr)
        print("Parsing PMF data: node {}".format(node))
    return ts_frames


def get_trj_path(node, format):
    return os.path.join(REORDER_TRJ_OUTPUT_DIR,
                        'pmf',
                        '{}.{}'.format(node, format))


def run_cpptraj(parm, ts_frames, format, output):

    for node, frames in ts_frames.items():
        if len(frames) == 0:
            continue

        print("Extracting TS frames: node {}, {} frames".format(node,
                                                                len(frames)))
        with open(CPPTRAJ_IN_NAME, 'w') as f:
            f.write('parm {}\n'.format(parm))
            trj_path = get_trj_path(node, format)
            frames_str = ','.join(map(str, frames))
            f.write('trajin {}\n'.format(trj_path))
            f.write('trajout {} onlyframes {}\n'.format(trj_path + '.tmp',
                                                    frames_str))

        os.system("cpptraj -i {} >> get_ts_frames.out".format(CPPTRAJ_IN_NAME))

    print("Combining")
    with open(CPPTRAJ_IN_NAME, 'w') as f:
        f.write('parm {}\n'.format(parm))
        for node, frames in ts_frames.items():
            if len(frames) == 0:
                continue
            trj_path = get_trj_path(node, format) + '.tmp'
            f.write('trajin {}\n'.format(trj_path))
        f.write('trajout {}\n'.format(output))

    os.system("cpptraj -i {} >> get_ts_frames.out".format(CPPTRAJ_IN_NAME))
    os.system("rm {}/*tmp".format(os.path.join(REORDER_TRJ_OUTPUT_DIR, 'pmf')))


def main():
    parser = argparse.ArgumentParser(
        description=""
        "Extracts frames close to TS from pmf trajectories obtained with "
        "reorder_trj.py.\n"
        "Must be executed from the working directory (the one where the STRING \n"
        "file is located)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("prmtop", help="Amber topology file")
    parser.add_argument("--pmf", help="PMF final to calculate the TS position")
    parser.add_argument("--format", help="Format of the pmf trajectories",
                        default='nc')
    parser.add_argument("--output", help="Name of the output trajectory",
                        default='ts.nc')
    parser.add_argument("--ts", type=float,
                        help="RC value to be used as TS", default=None)
    parser.add_argument("--thr", type=float,
                        help="TS threshold", default=0.05)
    args = parser.parse_args()

    ts = args.ts or get_ts(args.pmf)

    ts_frames = parse_pmf_dat(get_results_dir(), get_ntwx(), ts, args.thr)

    run_cpptraj(args.prmtop, ts_frames, args.format, args.output)


if __name__ == '__main__':
    main()
