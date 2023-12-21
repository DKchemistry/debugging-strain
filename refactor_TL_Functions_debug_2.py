import xml.etree.ElementTree as ET
from rdkit import Chem
import os
import numpy as np
from math import sqrt, atan2, pi
from math import ceil

tree = ET.parse("TL_2.1_VERSION_6.xml")
root = tree.getroot()


class TP_list(object):
    def __init__(self, indeces, angles, smarts, hc, methods, E, CI_l, CI_u, flags):
        self.indeces = indeces
        self.angles = angles
        self.smarts = smarts
        self.hc = hc
        self.methods = methods
        self.E = E
        self.CI_l = CI_l
        self.CI_u = CI_u
        self.flags = flags
        self.TP_indeces = [j for j in range(len(indeces))]

    def get_TPs(self, inds=None):
        if inds == None:
            inds = [j for j in range(len(self.indeces))]
        tps = []
        for j in inds:
            tps.append(
                [
                    self.TP_indeces[j],
                    self.E[j],
                    self.CI_l[j],
                    self.CI_u[j],
                    self.indeces[j],
                    self.angles[j],
                    self.smarts[j],
                    self.hc[j],
                    self.methods[j],
                    self.flags[j],
                ]
            )
        return tps

    def sum(self, cutoff=None):
        flagged = sum(self.flags)
        if flagged > 0:
            return [-1 * flagged, 0, 0]
        else:
            if cutoff == None:
                cutoff == 0
            ret = [0] * 3
            for i in range(len(self.E)):
                if self.E[i] >= cutoff:
                    ret[0] += self.E[i]
                    ret[1] += self.CI_l[i]
                    ret[2] += self.CI_u[i]
            return ret


def Mol2MolSupplier(file=None):
    names = []
    mols = {}
    with open(file, "r") as f:
        fileend = os.fstat(f.fileno()).st_size
        count = 0
        line = f.readline()
        while not f.tell() == fileend:
            if line.startswith("#") or line == "\n":
                line = f.readline()
            if line.startswith("@<TRIPOS>MOLECULE"):
                count += 1
                mol = []
                mol.append(line)
                line = f.readline()
                if line != "\n" and line.split()[0].strip() not in names:
                    name = line.split()[0].strip()
                else:
                    name = "mol2Number" + str(count)
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == fileend:
                        mol.append(line)
                        break
                block = ",".join(mol).replace(",", "")
                m = Chem.rdmolfiles.MolFromMol2Block(
                    block, sanitize=False, removeHs=False
                )
                names.append(name)
                mols[name] = m
    return (names, mols)


def custom_sdmolsupplier_H(file_path):
    if not file_path:
        raise ValueError("File path is required.")
    suppl_H = Chem.SDMolSupplier(file_path, removeHs=False)
    names = []
    mol = {}
    for mol_H in suppl_H:
        if mol_H is not None:
            name_H = (
                mol_H.GetProp("_Name")
                if mol_H.HasProp("_Name")
                else f"UnnamedMol_{len(names)}"
            )
            name_H = name_H.split()[0]
            names.append(name_H)
            mol[name_H] = mol_H
    return names, mol


def unit(a):
    return a / sqrt(np.dot(a, a))


def dihedral(a_1, a_2, a_3, a_4):
    b_1 = a_2 - a_1
    b_2 = a_3 - a_2
    b_3 = a_4 - a_3
    n_1 = unit(np.cross(b_1, b_2))
    n_2 = unit(np.cross(b_2, b_3))
    m = unit(np.cross(n_1, b_2))
    x = np.dot(n_1, n_2)
    y = np.dot(m, n_2)
    return -atan2(y, x) * 180 / pi


def ang_diff(theta_1, theta_2):
    if theta_1 < 0:
        theta_1 += 360
    if theta_2 < 0:
        theta_2 += 360
    del_theta = (theta_1 - theta_2) % 360
    if del_theta > 180:
        del_theta -= 360
    return del_theta


def tp_match(tp, hc, j, mol, pos, bi):
    smarts = tp.get("smarts")
    hist_E = []
    hist_l = []
    hist_u = []
    if tp.get("method") == "exact":
        for bin in tp.find("histogram_converted").findall("bin"):
            hist_E.append(float(bin.get("energy")))
            hist_l.append(float(bin.get("lower")))
            hist_u.append(float(bin.get("upper")))
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    for match in matches:
        if len(match) > 4:
            continue
        if (
            mol.GetAtomWithIdx(match[0]).GetSymbol() == "H"
            or mol.GetAtomWithIdx(match[3]).GetSymbol() == "H"
        ):
            continue
        r_1 = np.array(pos[match[0]])
        r_2 = np.array(pos[match[1]])
        r_3 = np.array(pos[match[2]])
        r_4 = np.array(pos[match[3]])
        theta = dihedral(r_1, r_2, r_3, r_4)
        if tp.get("method") == "exact":
            bin_num = ceil(theta / 10) + 17
            energy = (hist_E[bin_num] - hist_E[(bin_num + 35) % 36]) / 10.0 * (
                theta - (bin_num - 17) * 10
            ) + hist_E[bin_num]
            lower = (hist_l[bin_num] - hist_l[(bin_num + 35) % 36]) / 10.0 * (
                theta - (bin_num - 17) * 10
            ) + hist_l[bin_num]
            upper = (hist_u[bin_num] - hist_u[(bin_num + 35) % 36]) / 10.0 * (
                theta - (bin_num - 17) * 10
            ) + hist_u[bin_num]
            bi.append(
                [
                    list(match),
                    theta,
                    smarts,
                    hc,
                    "exact",
                    energy,
                    lower,
                    upper,
                    False,
                    j,
                ]
            )
        else:
            not_observed = True
            energy = 100
            for angle in tp.find("angleList").findall("angle"):
                theta_0 = float(angle.get("theta_0"))
                delta = ang_diff(theta, theta_0)
                if abs(delta) <= float(angle.get("tolerance2")):
                    beta_1 = float(angle.get("beta_1"))
                    beta_2 = float(angle.get("beta_2"))
                    energy = beta_1 * (delta**2) + beta_2 * (delta**4)
                    not_observed = False
                    break
            bi.append(
                [
                    list(match),
                    theta,
                    smarts,
                    hc,
                    "approximate",
                    energy,
                    energy,
                    energy,
                    not_observed,
                    j,
                ]
            )


def TL_lookup(mol):
    positions = mol.GetConformer().GetPositions()
    bond_info = []
    i = 0
    for HC in root.findall("hierarchyClass"):
        if HC.get("name") != "GG":
            for TP in HC.iter("torsionRule"):
                tp_match(TP, "specific", i, mol, positions, bond_info)
                i += 1
    for TP in root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
        tp_match(TP, "general", i, mol, positions, bond_info)
        i += 1
    for bond in bond_info:
        if bond[0][1] > bond[0][2]:
            bond[0].reverse()
            bond.append(True)
        else:
            bond.append(False)
    bond_info_red = [bond_info[0]]
    for j in range(1, len(bond_info)):
        atom_0 = bond_info[j][0][0]
        atom_1 = bond_info[j][0][1]
        atom_2 = bond_info[j][0][2]
        atom_3 = bond_info[j][0][3]
        unmatched = True
        for k in range(len(bond_info_red)):
            if (
                bond_info_red[k][0][0] == atom_0
                and bond_info_red[k][0][1] == atom_1
                and bond_info_red[k][0][2] == atom_2
                and bond_info_red[k][0][3] == atom_3
            ):
                unmatched = False
                if bond_info[j][9] < bond_info_red[k][9]:
                    bond_info_red[k] = bond_info[j]
                    break
        if unmatched:
            bond_info_red.append(bond_info[j])
    b_i_r = [bond_info_red[0]]
    for j in range(1, len(bond_info_red)):
        atom_1 = bond_info_red[j][0][1]
        atom_2 = bond_info_red[j][0][2]
        unmatched = True
        for k in range(len(b_i_r)):
            if b_i_r[k][0][1] == atom_1 and b_i_r[k][0][2] == atom_2:
                unmatched = False
                if bond_info_red[j][3][0] > b_i_r[k][3][0] or (
                    bond_info_red[j][5] > b_i_r[k][5]
                    and bond_info_red[j][3][0] == b_i_r[k][3][0]
                ):
                    b_i_r[k] = bond_info_red[j]
                    break
        if unmatched:
            b_i_r.append(bond_info_red[j])
    for bond in b_i_r:
        if bond[10]:
            bond[0].reverse()
    return TP_list(
        [bond[0] for bond in b_i_r],
        [bond[1] for bond in b_i_r],
        [bond[2] for bond in b_i_r],
        [bond[3] for bond in b_i_r],
        [bond[4] for bond in b_i_r],
        [bond[5] for bond in b_i_r],
        [bond[6] for bond in b_i_r],
        [bond[7] for bond in b_i_r],
        [bond[8] for bond in b_i_r],
    )
