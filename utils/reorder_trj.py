#!/usr/bin/env python

import argparse
import os
import shutil
import re
from operator import itemgetter
from itertools import groupby

TRAJECTORY_EXTENSION = 'nc'
SANDER_INPUT_EXTENSION = 'in'
OUTPUT_DIR = 'reorder_trj_results'
TMP_DIR = '.reorder_trj_tmp'
CPPTRAJ_IN_NAME = '.reorder_cpptraj_tmp.in'


def new_fragment(trj, node, start, stop):
    return {'trj': trj, 'node': node, 'start': start, 'stop': stop}


def step_to_frame(step, steps_per_frame, offset):
    return (step + offset) // steps_per_frame


def data_to_fragments(data, steps_per_frame, offset):
    ntrj = len(data[0][1:])
    fragments = []

    first_frame = step_to_frame(0, steps_per_frame, offset) + 1

    for trj in range(1, ntrj + 1):
        trj_fragments = []
        last_frame = -1
        for row in data:

            frame = step_to_frame(row[0], steps_per_frame, offset)
            if frame == last_frame or frame < first_frame:
                continue

            node = row[trj]
            if len(trj_fragments) == 0:
                trj_fragments.append(new_fragment(trj, node, first_frame, frame))
            elif trj_fragments[-1]['node'] == node:
                trj_fragments[-1]['stop'] = frame
            else:
                trj_fragments.append(new_fragment(trj, node, frame, frame))
            last_frame = frame

        fragments.append(trj_fragments)

    return fragments


def group_fragments_by_node(fragments):
    flat_fragments = [_ for trj_fragments in fragments for _ in trj_fragments]
    sorted_fragments = sorted(flat_fragments, key=itemgetter('node', 'start'))
    return groupby(sorted_fragments, key=itemgetter('node'))


def assign_fragments_filenames(fragments, path):
    for node, node_fragments in group_fragments_by_node(fragments):
        for index, fragment in enumerate(node_fragments):
            fragment['filename'] = os.path.join(path, '{}_{}.nc'.format(node, index))


def break_trj(in_filename, parm, trj_filename, fragments):
    with open(in_filename, 'w') as f:
        f.write('parm {}\n'.format(parm))
        f.write('trajin {}\n'.format(trj_filename))
        for fragment in fragments:
            f.write('trajout {filename} start {start} stop {stop}\n'
                    .format(**fragment))
    os.system('cpptraj -i {} >> break.tmp'.format(in_filename))


def combine_node_trj(in_filename, parm, node_trj_filename, node_fragments):
    with open(in_filename, 'w') as f:
        f.write('parm {}\n'.format(parm))
        for fragment in node_fragments:
            f.write('trajin {filename}\n'.format(**fragment))
        f.write('trajout {}\n'.format(node_trj_filename))
    os.system('cpptraj -i {} >> combine.tmp'.format(in_filename))


def mkdir(name):
    try:
        os.mkdir(name)
    except OSError:
        print('{} exists and will be removed'.format(name))
        shutil.rmtree(name)
        os.mkdir(name)


def move_dir_from_tmp(name):
    tmp_path = os.path.join(TMP_DIR, name)
    output_path = os.path.join(OUTPUT_DIR, name)
    os.system('mv {} {}'.format(tmp_path, output_path))


def reorder_trajectories(topology, filenames, trj_period,
                         offset, data, subdir, output_format):

    fragments = data_to_fragments(data, trj_period, offset)
    assign_fragments_filenames(fragments, TMP_DIR)

    for i, trj_fragments in enumerate(fragments):
        break_trj(CPPTRAJ_IN_NAME, topology, filenames[i], trj_fragments)
        print('Splitting {} of {} done'.format(i+1, len(filenames)))

    path = os.path.join(TMP_DIR, subdir)
    mkdir(path)

    for node, node_fragments in group_fragments_by_node(fragments):
        node_trj_path = os.path.join(path, str(node) + '.' + output_format)
        combine_node_trj(CPPTRAJ_IN_NAME, topology, node_trj_path, node_fragments)
        print('Combining {} of {} done'.format(node, len(filenames)))
    move_dir_from_tmp(subdir)
    if subdir == 'string':
        os.system("rm {}/*".format(TMP_DIR))


def read_data(filename):
    with open(filename) as f:
        return [[int(x) for x in line.strip().split()] for line in f]


def get_nnodes():
    """
    Assumes that for each node there is an 'in' file (1.in, 2.in etc) and
    there are no other 'in' files with numeric names.
    """
    regex = re.compile(r"([0-9]*)\."+SANDER_INPUT_EXTENSION)
    paths = filter(lambda x: regex.match(x) and os.path.isfile(x),
                   os.listdir('.'))
    return len(list(paths))


def trajectory_filenames(nnodes):
    return ['{}.{}'.format(i, TRAJECTORY_EXTENSION) for i in range(1, nnodes+1)]


def assert_file_exist(filename):
    if not os.path.exists(filename):
        print('{} not found. Exiting.'.format(filename))
        exit()


def parse_string():
    preparation_steps = None
    results_dir = '.'
    only_pmf = False
    assert_file_exist('STRING')
    with open('STRING') as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith('dir'):
                results_dir = line.split()[-1].strip('"')
            if line.startswith('preparation_steps'):
                preparation_steps = int(line.split()[-1])
            if line.startswith('only_PMF'):
                only_pmf = line.split()[-1].lower() == '.true.'
    return preparation_steps, results_dir, only_pmf


def parse_in(in_file):
    ntwx = None
    dt = None
    assert_file_exist(in_file)
    with open(in_file) as f:
        ntwx_regex = re.compile(r".*(^|[, \t]+)ntwx[ \t]*=[ \t]*([0-9]*)[\D]*")
        dt_regex = re.compile(r".*(^|[, \t]+)dt[ \t]*=[ \t]*([0-9\.]*)[\D]*")
        for line in f:
            ntwx_match = ntwx_regex.match(line.strip())
            if ntwx_match:
                ntwx = int(ntwx_match.groups()[1])
            dt_match = dt_regex.match(line.strip())
            if dt_match:
                dt = float(dt_match.groups()[1])
    return ntwx, dt


def parse_stop_string(stop_string):
    try:
        with open(stop_string) as f:
            return int(f.readlines()[-1].strip())
    except IOError:
        return None


def main():

    parser = argparse.ArgumentParser(
        description=""
        "Reorders trajectories generated by a string method calculation.\n"
        "Must be executed from the working directory (the one where the STRING \n"
        "file is located)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("prmtop", help="Amber topology file")
    parser.add_argument("--format", help="Format of the output trajectories",
                        default='nc')
    parser.add_argument("--nodes", type=int,
                        help="Number of nodes (use if automatic detection fails)")
    args = parser.parse_args()
    topology = args.prmtop

    # Parse STRING file to extract preparation_steps (default 2000),
    # results directory and only_PMF
    preparation_steps, results_dir, only_pmf = parse_string()

    # Parse first input file to know how often the frames were saved to the
    # trajectories (ntwx) and get the timestep (dt) to calculate
    # preparation_steps
    in_file = '1.{}'.format(SANDER_INPUT_EXTENSION)
    ntwx, dt = parse_in(in_file)
    if not ntwx:
        print('ntwx not set in {}, so no trajectories were written. Exiting.'
              .format(in_file))
        exit()
    if not preparation_steps:
        preparation_steps = int(1. / dt)

    # Parse convergence.dat (if found) to know after how many steps
    # the job switched to pmf
    try:
        with open(os.path.join(results_dir, 'convergence.dat')) as f:
            last_string_step = int(f.readlines()[-1].split()[0])
    except IOError:
        last_string_step = None

    has_string_data = os.path.exists(os.path.join(results_dir, 'plot.REX'))
    has_pmf_data = os.path.exists(os.path.join(results_dir, 'plot_final.REX'))
    if not (has_string_data or has_pmf_data):
        print('plot.REX or plot_final.REX not found in {}. Exiting.'
              .format(results_dir))
        exit()

    string_offset = preparation_steps
    pmf_offset = string_offset + last_string_step or 0
    if only_pmf:
        has_string = False
        has_pmf = True
    elif last_string_step:
        has_string = True
        has_pmf = True
    else:
        has_string = True
        has_pmf = False

    nnodes = args.nodes or get_nnodes()
    print("{} nodes".format(nnodes))

    trajectory_names = trajectory_filenames(nnodes)

    mkdir(TMP_DIR)
    mkdir(OUTPUT_DIR)

    if has_string:
        string_data = read_data(os.path.join(results_dir, 'plot.REX'))

        print("Extracting string trajectories...")
        reorder_trajectories(topology, trajectory_names, ntwx,
                             string_offset, string_data, 'string', args.format)

    if has_pmf:
        pmf_data = read_data(os.path.join(results_dir, 'plot_final.REX'))

        print("Extracting PMF trajectories...")
        reorder_trajectories(topology, trajectory_names, ntwx,
                             pmf_offset, pmf_data, 'pmf', args.format)

    shutil.rmtree(TMP_DIR)


if __name__ == '__main__':
    main()
