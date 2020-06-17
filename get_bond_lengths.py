#!/usr/bin/env python3

import os, sys

import openbabel
import pandas as pd


def read_mols(path):
	mols       = []
	mol_index  = []

	# Initialize the OpenBabel object that we will use later.
	obConversion = openbabel.OBConversion()
	obConversion.SetInFormat("xyz")

	for f in os.scandir(path):
		# Initialize an OBMol object
		mol = openbabel.OBMol()
		read_ok = obConversion.ReadFile(mol, f.path)
		if not read_ok:
			# There was an error reading the file
			raise Exception(f'Could not read file {f.path}')

		mol_name = f.name[:-4]
		mol_index.extend([mol_name] * mol.NumBonds()) # We need one entry per bond
		mols.append(mol)

	return mols, mol_index

def extract_data(mols, mol_index):
	bond_atom_0 = []
	bond_atom_1 = []
	bond_order  = []
	bond_length = []

	for mol in mols:
		# Extract bond information
		mol_bonds = openbabel.OBMolBondIter(mol) # get over all the bonds in the molecule
		for bond in mol_bonds:
			bond_atom_0.append(bond.GetBeginAtomIdx() - 1) # Must be 0-indexed
			bond_atom_1.append(bond.GetEndAtomIdx() - 1)
			bond_length.append(bond.GetLength())
			bond_order.append(bond.GetBondOrder())

	return bond_atom_0, bond_atom_1, bond_length, bond_order

if __name__ == '__main__':
	# path to xyz files
	path = sys.argv[1]

	# read in molecules
	mols, mol_index = read_mols(path)

	# extract data
	bond_atom_0, bond_atom_1, bond_length, bond_order = extract_data(mols, mol_index)

	# Put everything into a dataframe
	df = pd.DataFrame({'molecule_name': mol_index, 'atom_0': bond_atom_0, 'atom_1': bond_atom_1, 'order': bond_order, 'length': bond_length})

	df = df.sort_values(['molecule_name', 'atom_0', 'atom_1']).reset_index(drop=True)
	print(df)
