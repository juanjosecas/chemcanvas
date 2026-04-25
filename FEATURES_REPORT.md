# ChemCanvas – Additional Features Report

Analysis of the ChemCanvas codebase against the feature sets of
**ChemDraw**, **ACD/Labs**, and **Marvin Sketch** to identify
high-value, easily implementable additions.

All of the items below can be added with minimal disruption to the
existing architecture: they follow the same patterns already used by
`generateSmiles`, `readSmiles`, and `calculateProperties`.

---

## Category 1 – Molecular Properties & Analysis

| Feature | Basis | Difficulty |
|---------|-------|------------|
| **Molecular formula on canvas** – auto-render Hill formula below any selected molecule | atom symbol/charge already available; rdkit_properties already added | Low |
| **InChI / InChIKey generation** (ChemDraw, ACD) | RDKit `Chem.inchi.MolToInchi` / `InchiToInchiKey` | Low |
| **Molecular weight display in status bar** | RDKit already integrated | Low |
| **IUPAC name lookup** (Marvin, ChemDraw) | PubChem REST API via name → CID → iupac_name endpoint | Low |
| **pKa prediction** (ACD/pKa, Marvin) | `chemaxon` (if licensed) or `rdkit` + DMTA models; or call external web API | Medium |
| **Lipinski Ro5 / drug-likeness summary** | Pure RDKit calculation | Low |
| **Solubility estimate (ESOL)** | RDKit + ESOL equation (Delaney 2004) | Low |
| **Bioavailability radar** (ChemDraw 20+) | Requires matplotlib or Qt canvas drawing | Medium |

---

## Category 2 – Structure Operations

| Feature | Basis | Difficulty |
|---------|-------|------------|
| **Canonical SMILES normalisation** | RDKit `Chem.MolToSmiles(mol, canonical=True)` | Low |
| **Structure cleanup / 2-D coordinate regeneration** | RDKit `AllChem.Compute2DCoords` + feed back into existing `coords_generator` | Medium |
| **Tautomer enumeration** | RDKit `rdkit.Chem.MolStandardize.rdMolStandardize` | Medium |
| **Salt stripping / largest fragment** | RDKit `SaltRemover`, `LargestFragmentChooser` | Low |
| **Exact / substructure search in open molecules** | RDKit `HasSubstructMatch` + paper selection | Medium |
| **R-group decomposition** (ChemDraw) | RDKit `rdRGroupDecomposition` | High |
| **Reaction SMILES / RXN read–write** | ChemDraw, Marvin support; RDKit `AllChem` + `fileformat_rxn.py` following existing pattern | Medium |
| **SDF / MDL read-write** | `fileformat_molfile.py` already exists; extend to V3000 | Low |
| **CML (Chemical Markup Language) read-write** | Marvin native; `lxml` parsing, similar to `fileformat_mrv.py` | Low |

---

## Category 3 – Drawing & Rendering

| Feature | Basis | Difficulty |
|---------|-------|------------|
| **Atom numbering overlay** | Add numbered labels as `Text` objects following `calculateProperties` pattern | Low |
| **Bond length / angle labels** | Geometry already computed; render as `Text` near bonds | Low |
| **3-D conformation viewer** (Marvin 3D) | RDKit `AllChem.EmbedMolecule` + convert to z-coords; existing `transform_3D` support | High |
| **Color by atom type** (CPK colours) | Map element → colour; `atom.color` already supported | Low |
| **Highlight substructure** (ChemDraw) | Use existing `atom.set_selected` + `bond.set_selected` | Low |
| **Polymer/repeat-unit notation** (ChemDraw) | Extend `bracket.py` with polymer attributes | Medium |
| **Variable attachment point brackets** | Extend `bracket.py` | Medium |
| **Orbital/lone pair rendering toggle** | `atom.lonepairs` already stored; rendering already partially done | Low |

---

## Category 4 – Export & Integration

| Feature | Basis | Difficulty |
|---------|-------|------------|
| **Copy as SMILES to clipboard** | `QApplication.clipboard().setText(smiles)` one-liner | Low |
| **Copy as InChI to clipboard** | Same as above with InChI | Low |
| **Copy structure as bitmap/SVG** | `App.paper.getImage()` / `getSvg()` already exist; add Edit menu entry | Low |
| **Export to SDF with properties** | Extend `fileformat_molfile.py` with SD-tag writing via RDKit | Low |
| **Batch property calculation for all molecules** | Loop over `App.paper.objects` filtering by `class_name=="Molecule"` | Low |
| **PubChem / ChemSpider name search** (ChemDraw) | `template_manager.py` already queries PubChem; extend with name→SMILES | Low |
| **IUPAC → structure (name-to-structure)** | PubChem REST `compound/name/<name>/JSON` endpoint | Low |
| **CAS number lookup** | PubChem REST | Low |

---

## Category 5 – User Interface

| Feature | Basis | Difficulty |
|---------|-------|------------|
| **Properties panel (dockable)** | `QDockWidget` on the right side; update on selection change | Medium |
| **Structure validation warnings** | RDKit `Chem.SanitizeMol` error reporting | Low |
| **Undo history viewer** | `undo_manager.py` already records names; display in a `QListWidget` | Low |
| **Dark-mode aware chemical colours** | `App.dark_mode` flag already exists; adjust bond/atom default colours | Low |
| **Keyboard shortcut for tools** | `QAction.setShortcut()` in `tools_template`; several tools lack shortcuts | Low |

---

## Implementation Notes

- Every item in categories 1–2 that uses RDKit can follow the same pattern
  as `calculateProperties` in `main.py` + `rdkit_properties.py`:
  pure function → method in `Window` → action in `menuTools`.
- Items requiring new drawing primitives should extend `drawing_parents.py`
  and follow the `DrawableObject` contract.
- Network-based features (PubChem queries) should use the `QThread`/signal
  pattern already demonstrated by the `UpdateChecker` class in `widgets.py`.
