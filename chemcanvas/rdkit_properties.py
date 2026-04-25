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


# ---------------------------------------------------------------------------
# InChI / InChIKey
# ---------------------------------------------------------------------------

def generate_inchi(smiles):
    """Generate InChI and InChIKey from a SMILES string.

    Returns a dict with keys ``inchi`` and ``inchikey`` on success, or
    ``None`` if RDKit is not available, SMILES is invalid, or InChI
    generation fails.
    """
    if not _rdkit_available:
        return None
    try:
        from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
    except ImportError:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    inchi = MolToInchi(mol)
    if inchi is None:
        return None
    inchikey = InchiToInchiKey(inchi) or ""
    return {"inchi": inchi, "inchikey": inchikey}


# ---------------------------------------------------------------------------
# Canonical SMILES
# ---------------------------------------------------------------------------

def canonical_smiles(smiles):
    """Return the RDKit canonical SMILES for the given SMILES string.

    Returns ``None`` if RDKit is not available or the SMILES is invalid.
    """
    if not _rdkit_available:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


# ---------------------------------------------------------------------------
# Lipinski Ro5
# ---------------------------------------------------------------------------

def lipinski_ro5(smiles):
    """Check Lipinski's Rule of Five for a SMILES string.

    Returns a dict with keys:
        mol_weight, logp, hbd, hba, violations, drug_like

    ``violations`` is the number of rules broken (drug-like if ≤ 1).
    Returns ``None`` if RDKit is not available or the SMILES is invalid.
    """
    if not _rdkit_available:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mw   = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return {
        "mol_weight": mw,
        "logp":       logp,
        "hbd":        hbd,
        "hba":        hba,
        "violations": violations,
        "drug_like":  violations <= 1,
    }


def format_lipinski_html(props):
    """Format Lipinski Ro5 results as an HTML string."""
    if props is None:
        return ""
    verdict = ("<font color='green'><b>Drug-like</b></font>"
               if props["drug_like"] else
               "<font color='red'><b>Not drug-like</b></font>")
    lines = [
        verdict,
        "<b>MW:</b> {mol_weight} g/mol  (≤ 500)",
        "<b>LogP:</b> {logp}  (≤ 5)",
        "<b>HBD:</b> {hbd}  (≤ 5)",
        "<b>HBA:</b> {hba}  (≤ 10)",
        "<b>Violations:</b> {violations} / 4",
    ]
    return "<br>".join(line.format(**props) for line in lines)


def format_lipinski_plain(props):
    """Format Lipinski Ro5 results as plain text."""
    if props is None:
        return ""
    verdict = "Drug-like (Ro5)" if props["drug_like"] else "NOT drug-like (Ro5 violated)"
    lines = [
        verdict,
        "MW:         {mol_weight} g/mol  (rule: <= 500)",
        "LogP:       {logp}  (rule: <= 5)",
        "HBD:        {hbd}  (rule: <= 5)",
        "HBA:        {hba}  (rule: <= 10)",
        "Violations: {violations} / 4",
    ]
    return "\n".join(line.format(**props) for line in lines)


# ---------------------------------------------------------------------------
# Salt stripping / largest fragment
# ---------------------------------------------------------------------------

def strip_salts(smiles):
    """Return the SMILES of the largest organic fragment (salt stripping).

    Returns the canonical SMILES of the largest fragment, or ``None`` if
    RDKit is not available or the SMILES is invalid.
    """
    if not _rdkit_available:
        return None
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    chooser = rdMolStandardize.LargestFragmentChooser()
    fragment = chooser.choose(mol)
    return Chem.MolToSmiles(fragment, canonical=True)


# ---------------------------------------------------------------------------
# Structure validation
# ---------------------------------------------------------------------------

def validate_structure(smiles):
    """Validate a SMILES string using RDKit sanitisation.

    Returns a list of warning strings.  An empty list means the structure
    is valid.  Returns ``None`` if RDKit is not available.
    """
    if not _rdkit_available:
        return None

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return ["Could not parse SMILES string."]

    errors = []
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        errors.append(str(e))
    return errors
