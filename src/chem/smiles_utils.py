"""RDKit-backed SMILES generation from ChemCanvas-style molecular graphs."""

from rdkit import Chem
from rdkit.Chem import Draw


_BOND_TYPE_MAP = {
    "single": Chem.BondType.SINGLE,
    "wedge": Chem.BondType.SINGLE,
    "hashed_wedge": Chem.BondType.SINGLE,
    "wavy": Chem.BondType.SINGLE,
    "bold": Chem.BondType.SINGLE,
    "double": Chem.BondType.DOUBLE,
    "E_or_Z": Chem.BondType.DOUBLE,
    "bold2": Chem.BondType.DOUBLE,
    "triple": Chem.BondType.TRIPLE,
    "delocalized": Chem.BondType.AROMATIC,
    "coordinate": Chem.BondType.DATIVE,
}


class _TestNamespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def _graph_atoms(graph):
    atoms = getattr(graph, "atoms", None)
    if atoms is None:
        atoms = getattr(graph, "vertices", None)
    if atoms is None:
        raise ValueError("Graph has no atoms or vertices")
    return list(atoms)


def _graph_bonds(graph):
    bonds = getattr(graph, "bonds", None)
    if bonds is None:
        bonds = getattr(graph, "edges", None)
    if bonds is None:
        raise ValueError("Graph has no bonds or edges")
    return list(bonds)


def _atom_symbol(atom):
    symbol = getattr(atom, "symbol", None)
    if not symbol:
        atomic_num = getattr(atom, "atomic_num", None)
        if atomic_num is None:
            atomic_num = getattr(atom, "atomic_number", None)
        if atomic_num:
            symbol = Chem.GetPeriodicTable().GetElementSymbol(int(atomic_num))
    if not symbol:
        raise ValueError("Atom is missing an element symbol or atomic number")
    return symbol


def _formal_charge(atom):
    charge = getattr(atom, "charge", 0)
    properties = getattr(atom, "properties_", {}) or {}
    if not charge and "charge" in properties:
        charge = properties["charge"]
    return int(charge or 0)


def _bond_atoms(bond):
    atoms = getattr(bond, "atoms", None)
    if atoms is None:
        atoms = getattr(bond, "vertices", None)
    if not atoms or len(atoms) != 2:
        raise ValueError("Bond does not connect exactly two atoms")
    return atoms[0], atoms[1]


def _bond_order(bond):
    order = getattr(bond, "order", None)
    if order is None:
        bond_type = getattr(bond, "type", "single")
        return {"single": 1, "double": 2, "triple": 3, "delocalized": 1.5}.get(bond_type)
    return order


def _bond_type_name(bond):
    return getattr(bond, "type", None) or {1: "single", 2: "double", 3: "triple", 1.5: "delocalized"}.get(_bond_order(bond), "single")


def _is_atom_aromatic(atom):
    properties = getattr(atom, "properties_", {}) or {}
    return bool(properties.get("aromatic") or getattr(atom, "aromatic", False) or getattr(atom, "is_aromatic", False))


def _is_bond_aromatic(bond):
    properties = getattr(bond, "properties_", {}) or {}
    return bool(
        _bond_type_name(bond) == "delocalized"
        or properties.get("aromatic")
        or properties.get("delocalized")
        or getattr(bond, "aromatic", False)
        or getattr(bond, "is_aromatic", False)
    )


def _find_bond_between(atom1, atom2):
    for bond in getattr(atom1, "bonds", getattr(atom1, "edges", [])):
        try:
            b_atom1, b_atom2 = _bond_atoms(bond)
        except ValueError:
            continue
        if {b_atom1, b_atom2} == {atom1, atom2}:
            return bond
    return None


def _canonical_cycle(cycle):
    rotations = []
    values = list(cycle)
    for seq in (values, list(reversed(values))):
        for i in range(len(seq)):
            rotations.append(tuple(seq[i:] + seq[:i]))
    return min(rotations, key=lambda atoms: tuple(id(atom) for atom in atoms))


def _six_membered_cycles(atoms):
    cycles = set()
    atom_set = set(atoms)
    for start in atoms:
        stack = [(start, [start])]
        while stack:
            current, path = stack.pop()
            if len(path) > 6:
                continue
            for neighbor in getattr(current, "neighbors", []):
                if neighbor not in atom_set:
                    continue
                if neighbor is start and len(path) == 6:
                    cycles.add(_canonical_cycle(path))
                elif neighbor not in path and len(path) < 6:
                    stack.append((neighbor, path + [neighbor]))
    return [list(cycle) for cycle in cycles]


def _is_alternating_six_ring(ring):
    bonds = []
    orders = []
    for i, atom in enumerate(ring):
        bond = _find_bond_between(atom, ring[(i + 1) % len(ring)])
        if bond is None:
            return None
        bonds.append(bond)
        orders.append(_bond_order(bond))
    if sorted(orders) != [1, 1, 1, 2, 2, 2]:
        return None
    for i in range(6):
        if orders[i] == orders[(i + 1) % 6]:
            return None
    return bonds


def _aromatic_sets(atoms, bonds):
    aromatic_atoms = {atom for atom in atoms if _is_atom_aromatic(atom)}
    aromatic_bonds = {bond for bond in bonds if _is_bond_aromatic(bond)}

    for bond in aromatic_bonds:
        atom1, atom2 = _bond_atoms(bond)
        aromatic_atoms.update((atom1, atom2))

    for ring in _six_membered_cycles(atoms):
        ring_bonds = _is_alternating_six_ring(ring)
        if ring_bonds:
            aromatic_atoms.update(ring)
            aromatic_bonds.update(ring_bonds)

    return aromatic_atoms, aromatic_bonds


def graph_to_rdkit_mol(graph):
    """Convert a ChemCanvas-style graph to an RDKit Mol."""
    atoms = _graph_atoms(graph)
    bonds = _graph_bonds(graph)
    aromatic_atoms, aromatic_bonds = _aromatic_sets(atoms, bonds)
    atom_to_idx = {}
    rw_mol = Chem.RWMol()
    periodic_table = Chem.GetPeriodicTable()

    for atom in atoms:
        symbol = _atom_symbol(atom)
        atomic_num = periodic_table.GetAtomicNumber(symbol)
        if atomic_num <= 0:
            raise ValueError("Unsupported atom symbol: %s" % symbol)

        rd_atom = Chem.Atom(int(atomic_num))
        rd_atom.SetFormalCharge(_formal_charge(atom))
        isotope = getattr(atom, "isotope", None)
        if isotope:
            rd_atom.SetIsotope(int(isotope))
        if atom in aromatic_atoms:
            rd_atom.SetIsAromatic(True)
        if not getattr(atom, "auto_hydrogens", True):
            rd_atom.SetNumExplicitHs(int(getattr(atom, "hydrogens", 0) or 0))
            rd_atom.SetNoImplicit(True)
        atom_to_idx[atom] = rw_mol.AddAtom(rd_atom)

    for bond in bonds:
        atom1, atom2 = _bond_atoms(bond)
        if atom1 not in atom_to_idx or atom2 not in atom_to_idx:
            raise ValueError("Bond references an atom outside the graph")
        if bond in aromatic_bonds:
            rd_bond_type = Chem.BondType.AROMATIC
        else:
            bond_type_name = _bond_type_name(bond)
            if bond_type_name not in _BOND_TYPE_MAP:
                raise ValueError("Unsupported bond type for SMILES export: %s" % bond_type_name)
            rd_bond_type = _BOND_TYPE_MAP[bond_type_name]
        rw_mol.AddBond(atom_to_idx[atom1], atom_to_idx[atom2], rd_bond_type)
        rd_bond = rw_mol.GetBondBetweenAtoms(atom_to_idx[atom1], atom_to_idx[atom2])
        if rd_bond_type == Chem.BondType.AROMATIC:
            rd_bond.SetIsAromatic(True)

    return rw_mol.GetMol()


def sanitize_mol(mol):
    """Sanitize an RDKit molecule and return (mol, error)."""
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return mol, str(exc)
    return mol, None


def generate_valid_smiles(graph, canonical=True):
    try:
        mol = graph_to_rdkit_mol(graph)
        mol, error = sanitize_mol(mol)
        if error:
            return {"smiles": None, "valid": False, "error": error}
        smiles = Chem.MolToSmiles(mol, canonical=canonical)
        mol2 = Chem.MolFromSmiles(smiles)
        if mol2 is None:
            return {"smiles": smiles, "valid": False, "error": "RDKit could not parse generated SMILES"}
        return {"smiles": smiles, "valid": True, "error": None}
    except Exception as exc:
        return {"smiles": None, "valid": False, "error": str(exc)}


def draw_mol_debug(mol):
    return Draw.MolToImage(mol)


def _test_atom(symbol, charge=0, hydrogens=0, auto_hydrogens=True):
    atom = _TestNamespace(
        symbol=symbol,
        charge=charge,
        isotope=None,
        hydrogens=hydrogens,
        auto_hydrogens=auto_hydrogens,
        properties_={},
        neighbors=[],
        bonds=[],
    )
    return atom


def _test_bond(atom1, atom2, bond_type="single"):
    bond = _TestNamespace(type=bond_type, atoms=[atom1, atom2], vertices=[atom1, atom2], properties_={})
    atom1.neighbors.append(atom2)
    atom2.neighbors.append(atom1)
    atom1.bonds.append(bond)
    atom2.bonds.append(bond)
    return bond


def _test_graph(atoms, bonds):
    return _TestNamespace(atoms=atoms, vertices=atoms, bonds=set(bonds), edges=set(bonds))


def test_smiles_generation():
    methane = _test_graph([_test_atom("C")], [])
    result = generate_valid_smiles(methane)
    assert result["valid"] is True and result["error"] is None

    c1 = _test_atom("C")
    c2 = _test_atom("C")
    o = _test_atom("O")
    ethanol = _test_graph([c1, c2, o], [_test_bond(c1, c2), _test_bond(c2, o)])
    result = generate_valid_smiles(ethanol)
    assert result["valid"] is True and result["error"] is None

    ring_atoms = [_test_atom("C") for _ in range(6)]
    ring_bonds = []
    for i in range(6):
        ring_bonds.append(_test_bond(ring_atoms[i], ring_atoms[(i + 1) % 6], "double" if i % 2 == 0 else "single"))
    benzene = _test_graph(ring_atoms, ring_bonds)
    result = generate_valid_smiles(benzene)
    assert result["valid"] is True and result["error"] is None

    ammonium = _test_graph([_test_atom("N", charge=1, hydrogens=4, auto_hydrogens=False)], [])
    result = generate_valid_smiles(ammonium)
    assert result["valid"] is True and result["error"] is None
