# -*- coding: utf-8 -*-
# This file is a part of ChemCanvas Program which is GNU GPLv3 licensed
# Copyright (C) 2022-2026 Arindam Chaudhuri <arindamsoft94@gmail.com>
#
# RDKit physicochemical property calculations
# All functions are module-level (no classes required).

_rdkit_available = False

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    _rdkit_available = True
except ImportError:
    pass


def rdkit_available():
    """Return True if RDKit is installed and importable."""
    return _rdkit_available


def calculate_properties(smiles):
    """Calculate common physicochemical properties from a SMILES string.

    Returns a dict with the following keys on success:
        formula, mol_weight, exact_mass, logp, tpsa,
        hbd, hba, rotatable_bonds, rings, aromatic_rings

    Returns None if RDKit is not available or if the SMILES is invalid.
    """
    if not _rdkit_available:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "formula":         rdMolDescriptors.CalcMolFormula(mol),
        "mol_weight":      round(Descriptors.MolWt(mol), 4),
        "exact_mass":      round(Descriptors.ExactMolWt(mol), 6),
        "logp":            round(Descriptors.MolLogP(mol), 2),
        "tpsa":            round(Descriptors.TPSA(mol), 2),
        "hbd":             rdMolDescriptors.CalcNumHBD(mol),
        "hba":             rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings":           rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings":  rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def format_properties_html(props):
    """Format a property dict as an HTML string suitable for QGraphicsTextItem.

    Each property is shown on its own line using a small, readable font.
    """
    if props is None:
        return ""

    lines = [
        "<b>Formula:</b> {formula}",
        "<b>MW:</b> {mol_weight} g/mol",
        "<b>Exact Mass:</b> {exact_mass} g/mol",
        "<b>LogP:</b> {logp}",
        "<b>TPSA:</b> {tpsa} Å²",
        "<b>HBD / HBA:</b> {hbd} / {hba}",
        "<b>Rot. Bonds:</b> {rotatable_bonds}",
        "<b>Rings (arom.):</b> {rings} ({aromatic_rings})",
    ]
    return "<br>".join(line.format(**props) for line in lines)


def format_properties_plain(props):
    """Format a property dict as plain text (one property per line)."""
    if props is None:
        return ""

    lines = [
        "Formula:        {formula}",
        "MW:             {mol_weight} g/mol",
        "Exact Mass:     {exact_mass} g/mol",
        "LogP:           {logp}",
        "TPSA:           {tpsa} Å²",
        "HBD / HBA:      {hbd} / {hba}",
        "Rot. Bonds:     {rotatable_bonds}",
        "Rings (arom.):  {rings} ({aromatic_rings})",
    ]
    return "\n".join(line.format(**props) for line in lines)
