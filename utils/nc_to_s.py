#!/usr/bin/env python3

import argparse

import numpy as np
from netCDF4 import Dataset
import molmod as mm

# HPC specific, may not be needed
import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'


class CVs:
    CV_TYPE_DICT = {
        'BOND': lambda xyz: mm.bond_length(xyz)[0],
        'ANGLE': lambda xyz: mm.bend_angle(xyz)[0] / mm.deg,
        'DIHEDRAL': lambda xyz: mm.dihed_angle(xyz)[0] / mm.deg,
        'PPLANE': lambda xyz: CVs.eval_pplane(np.roll(xyz, -3))
    }

    def __init__(self, cvs_filename):
        with open(cvs_filename) as f:
            lines = f.readlines()
        self.cv_type_list = [self.extract_cv_type(_) for _ in lines
                             if _.strip().lower().startswith('colvar_type')]
        self.atoms_list = [self.extract_atoms(_) for _ in lines
                           if _.strip().lower().startswith('atoms')]

    def __call__(self, xyz):
        return np.array([self.get_cv(*_, xyz) for _ in
                         zip(self.cv_type_list, self.atoms_list)])

    @staticmethod
    def extract_cv_type(line):
        return line.split('=')[1].strip().strip('"').strip("'")

    @staticmethod
    def extract_atoms(line):
        return [int(_) - 1 for _ in line.split('=')[1].strip().split(',')]

    @classmethod
    def get_cv(cls, cv_type, atoms, xyz):
        atoms_xyz = xyz[atoms, :]
        return cls.CV_TYPE_DICT[cv_type](atoms_xyz)

    @staticmethod
    def eval_pplane(xyz):
        return mm.opbend_dist(xyz)[0] / mm.angstrom


class PathCV:

    def __init__(self, cvs_filename, def_filename):
        self.cvs_calculator = CVs(cvs_filename)
        with open(def_filename) as f:
            raw_ncvs, raw_npoints, raw_lamb = next(f).split()
            self.ncvs = int(raw_ncvs)
            self.npoints = int(raw_npoints)
            self.lamb = float(raw_lamb)
            self.arc = np.array([float(next(f))
                                 for _ in range(self.npoints)])
            self.theta = np.array([self.read_vector(f)
                                   for _ in range(self.npoints)])
            self.M = np.array([self.read_matrix(f, self.ncvs)
                               for _ in range(self.npoints)])
            self.dihedrals = [_ == 'DIHEDRAL'
                              for _ in self.cvs_calculator.cv_type_list]

    def __call__(self, xyz):
        cvs = self.cvs_calculator(xyz)
        dtheta = self.theta - cvs
        dtheta[:, self.dihedrals] = self.map_dihedral(dtheta[:, self.dihedrals])
        d = np.array([self.dist_M(*_) for _ in zip(dtheta, self.M)])
        w = np.exp(-d * self.lamb)
        return w @ self.arc / np.sum(w), cvs

    @classmethod
    def read_matrix(cls, f, n):
        return np.array([cls.read_vector(f) for _ in range(n)])

    @staticmethod
    def read_vector(f):
        return np.array(list(map(float, next(f).split())))

    @staticmethod
    def dist_M(r, M):
        return np.sqrt(r @ M @ r)

    @staticmethod
    def map_dihedral(dx):
        return (dx + 180) % 360 - 180


def nc_to_pathcv(nc, cvs, pathcv_def):
    f = Dataset(nc)
    xyz = f['coordinates']
    pathcv_calculator = PathCV(cvs, pathcv_def)
    result = list(map(np.array, zip(*[pathcv_calculator(_) for _ in xyz])))
    f.close()
    return result


parser = argparse.ArgumentParser(
    description="Calculate pathCV and CVs values from NetCDF trajectory"
)
parser.add_argument("cvs", help="CVs file")
parser.add_argument("trj", help="NetCDF trajectory file")
parser.add_argument("pathcv", help="PathCV definition")
parser.add_argument("out", help="Output file")
args = parser.parse_args()

s, cvs = nc_to_pathcv(args.trj, args.cvs, args.pathcv)
np.savetxt(args.out, np.hstack([s[:, None], cvs]), '%12.6f')
