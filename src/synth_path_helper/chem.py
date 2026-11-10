"""Low-level RDKit helpers: parsing input, canonical SMILES and product clean-up."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

# RDKit prints a parse error to stderr for every bad SMILES. We report errors
# ourselves, so silence the C++ logger.
RDLogger.DisableLog("rdApp.*")

DATA_DIR = Path(__file__).resolve().parent / "data"

MoleculeInput = Union[str, Chem.Mol]


class InvalidMoleculeError(ValueError):
    """Raised when an input cannot be turned into a valid RDKit molecule."""


@lru_cache(maxsize=1)
def common_names() -> Dict[str, str]:
    """Return the built-in name -> SMILES dictionary (keys are lower-case)."""
    with open(DATA_DIR / "common_compounds.json", encoding="utf-8") as fh:
        data = json.load(fh)
    return {k.lower(): v for k, v in data["compounds"].items()}


def _normalise(mol: Chem.Mol) -> Chem.Mol:
    """Round-trip through canonical SMILES.

    This drops explicit-H flags on atoms such as ``[CH3]`` that the user may
    have typed. Otherwise those flags would survive into reaction products
    and give wrong hydrogen counts.
    """
    smiles = Chem.MolToSmiles(mol)
    clean = Chem.MolFromSmiles(smiles)
    if clean is None:  # pragma: no cover - canonical SMILES should always parse
        raise InvalidMoleculeError(f"RDKit could not re-read its own SMILES {smiles!r}")
    return clean


def parse_molecule(value: MoleculeInput) -> Chem.Mol:
    """Parse a SMILES string, a common name, or an existing ``Chem.Mol``.

    Names are looked up first in the small built-in dictionary (e.g.
    ``"cyclohexanol"``); anything else is treated as SMILES.
    """
    return _normalise(_parse_raw(value))


def _parse_raw(value: MoleculeInput) -> Chem.Mol:
    if isinstance(value, Chem.Mol):
        return Chem.Mol(value)
    if not isinstance(value, str):
        raise InvalidMoleculeError(f"Expected a SMILES string or RDKit Mol, got {type(value).__name__}")
    text = value.strip()
    if not text:
        raise InvalidMoleculeError("Empty input")
    smiles = common_names().get(text.lower(), text)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidMoleculeError(
            f"Could not parse {value!r} as SMILES (and it is not a known compound name)"
        )
    if mol.GetNumAtoms() == 0:
        raise InvalidMoleculeError("Input contains no atoms")
    return mol


def parse_mixture(values: Union[MoleculeInput, Iterable[MoleculeInput]]) -> List[Chem.Mol]:
    """Parse one or more inputs and return the individual molecules (components).

    A dot-separated SMILES such as ``"CC(=O)O.CCO"`` gives two components, as
    does ``["CC(=O)O", "CCO"]``. Components are returned in the order typed.
    """
    if isinstance(values, (str, Chem.Mol)):
        values = [values]
    components: List[Chem.Mol] = []
    for value in values:
        mol = _parse_raw(value)
        for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True):
            components.append(_normalise(frag))
    if not components:
        raise InvalidMoleculeError("No molecules given")
    return components


def canonical_smiles(mol: Chem.Mol, isomeric: bool = True) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def canonicalize(smiles: str, isomeric: bool = True) -> str:
    """Canonicalise a SMILES string; raises :class:`InvalidMoleculeError`."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidMoleculeError(f"Invalid SMILES {smiles!r}")
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def formula(mol: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def mol_weight(mol: Chem.Mol) -> float:
    return round(Descriptors.MolWt(mol), 2)


def finalize_product(raw: Chem.Mol) -> Optional[Chem.Mol]:
    """Clean up a raw ``RunReactants`` product and sanitise it.

    RDKit copies mapped atoms from the reactant, including any 'no implicit
    H' flag. For atoms whose bonding changed in the reaction we let RDKit
    recompute the hydrogens and drop stale chirality. Returns ``None`` if the
    product is not chemically valid.

    The returned molecule keeps the reaction atom properties (``old_mapno``)
    that :mod:`synth_path_helper.selectivity` uses.
    """
    mol = Chem.Mol(raw)
    for atom in mol.GetAtoms():
        if not atom.HasProp("old_mapno"):
            continue
        if not atom.GetIsAromatic():
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
        if atom.HasProp("_ReactionDegreeChanged"):
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
    except Exception:  # RDKit raises several exception types for bad valences/aromaticity
        return None
    # Final validity check: the canonical SMILES has to parse again.
    if Chem.MolFromSmiles(Chem.MolToSmiles(mol)) is None:
        return None
    return mol
