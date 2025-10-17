"""
this package is used to generate CG_RNA model
Athour: Shanlong Li
Date: Nov 13, 2023
"""

import sys
import argparse
import numpy as np
from pathlib import Path

AAs = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE',
       'L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}

def printcg(atoms, file):
    for atom in atoms:
        file.write('{}  {:5d} {:>2}   {} {}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {}\n'.format(atom[0],int(atom[1]),atom[2],atom[3],atom[4],int(atom[5]),atom[6],atom[7],atom[8],atom[9],atom[10], atom[11]))

def read_map(seq):
    atoms = []
    filename = "map/protein/"+seq+".map"
    f_map = Path(__file__).parent / filename
    data = np.genfromtxt(f_map, dtype=None, names=('index', 'name', 'rx', 'ry', 'rz'), encoding='utf-8')
    for l in data:
        atom = ['ATOM', l[0], l[1], AAs[seq], 'X', 1, l[2], l[3], l[4], 1.00, 0.00, 'PRO']
        atoms.append(atom)
    return atoms

def transform(CA, CO, atoms, theta=120.0, legnth=1.92):
    nCA, nCO = np.array([atoms[0][6], atoms[0][7], atoms[0][8]]), np.array([atoms[1][6], atoms[1][7], atoms[1][8]])   # positions of CA and CO for next residue
    def norm(v):
        return v / np.linalg.norm(v)
    
    def rotation_matrix(axis, theta=theta):
        axis = norm(axis)
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        I = np.eye(3)
        return I + np.sin(theta) * K + (1 - np.cos(theta)) * np.dot(K, K)
    
    CACO = np.array(CO) - np.array(CA)
    CACO_norm = norm(CACO)
    new = CO + legnth * CACO_norm
    atoms[0][6], atoms[0][7], atoms[0][8] = new[0], new[1], new[2] # update CA position

    COCA = CO - new
    COCA_norm = norm(COCA)
    
    nCACO = nCO - nCA
    rot_axis = np.cross(COCA_norm, nCACO)
    rot_axis = norm(rot_axis)
    angle = np.radians(theta)
    R = rotation_matrix(rot_axis, angle)

    v_new = R @ COCA
    new_CO = new + norm(v_new) * np.linalg.norm(nCACO)

    v1, v2 = norm(nCACO), norm(new_CO - new)
    cross = np.cross(v1, v2)
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    if np.linalg.norm(cross) < 1e-8:
        R2 = np.eye(3)
    else:
        rot_axis2 = norm(cross)
        angle2 = np.arccos(dot)
        R2 = rotation_matrix(rot_axis2, angle2)

    for atom in atoms[1:]:
        vec = np.array([atom[6], atom[7], atom[8]]) - nCA
        new_atom = new + R2 @ (vec)
        atom[6] = new_atom[0]
        atom[7] = new_atom[1]
        atom[8] = new_atom[2]

    def rotate_around_axis(point, axis_point1, axis_point2, angle_rad):
        axis_dir = norm(axis_point2 - axis_point1)
        R = rotation_matrix(axis_dir, angle_rad)
        vec = point - axis_point1
        rotated = R @ vec
        return axis_point1 + rotated

    
    if len(atoms) > 2:
        CA_new = np.array([atoms[0][6], atoms[0][7], atoms[0][8]])
        CO_new = np.array([atoms[1][6], atoms[1][7], atoms[1][8]])
        CB_new = np.array([atoms[2][6], atoms[2][7], atoms[2][8]])
        if np.linalg.norm(CB_new - CO) <= 2.92:
            for angle_deg_step in np.linspace(10, 360, 36):  # 10° steps
                angle_step_rad = np.radians(angle_deg_step)
                F_rotated = rotate_around_axis(CB_new, CA_new, CO_new, angle_step_rad)
                if np.linalg.norm(F_rotated - CO) > 2.92:
                    for atom in atoms[2:]:
                        point = np.array([atom[6], atom[7], atom[8]])
                        rotated_point = rotate_around_axis(point, CA_new, CO_new, angle_step_rad)
                        atom[6], atom[7], atom[8] = rotated_point[0], rotated_point[1], rotated_point[2]
                    break

    return atoms

def build(seqs, out):
    with open(out, 'w') as f:
        print('REMARK  iCon Protein Model', file=f)
        print('REMARK  CREATE BY ProteinBuilder/SHANLONG LI', file=f)
        idx = 0
        res = 0
        for i, seq in enumerate(seqs):
            atoms = read_map(seq)
            if i == 0:
                atoms = atoms
            else:
                atoms = transform(CA, CO, atoms)
            CA = [atoms[0][6], atoms[0][7], atoms[0][8]]
            CO = [atoms[1][6], atoms[1][7], atoms[1][8]]
            
            for atom in atoms:
                atom[1] += idx
                atom[5] += res
            idx += len(atoms)
            res += 1
            printcg(atoms, f)
        print('END', file=f)
