"""Multi-step route search: breadth-first search over one-component transformations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from rdkit import Chem

from .chem import MoleculeInput, canonical_smiles, parse_molecule
from .functional_groups import detect_groups
from .predictor import Prediction, predict
from .transformations import Transformation, load_transformations


@dataclass(frozen=True)
class RouteStep:
    transformation: Transformation
    reactant: str
    product: str           # the fragment carried on to the next step
    full_product: str      # everything the step gives (e.g. both ozonolysis fragments)
    warnings: Tuple[str, ...] = ()
    co_reactants: Tuple[str, ...] = ()  # partner molecules used by two-component steps

    def to_dict(self) -> Dict[str, Any]:
        t = self.transformation
        return {
            "transformation": t.id,
            "name": t.name,
            "reactant": self.reactant,
            "co_reactants": list(self.co_reactants),
            "product": self.product,
            "full_product": self.full_product,
            "reagents": list(t.reagents),
            "catalysts": list(t.catalysts),
            "solvents": list(t.solvents),
            "warnings": list(self.warnings),
        }


@dataclass
class Route:
    start: str
    target: str
    steps: List[RouteStep] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.steps)

    @property
    def warning_count(self) -> int:
        return sum(len(s.warnings) for s in self.steps)

    def to_dict(self) -> Dict[str, Any]:
        return {"start": self.start, "target": self.target, "steps": [s.to_dict() for s in self.steps]}


def _key(smiles: str) -> str:
    """Stereo-insensitive canonical key used to compare molecules."""
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=False) if mol is not None else smiles


def _candidate_predictions(
    mol: Chem.Mol,
    groups: set,
    library: Sequence[Transformation],
    partners: Sequence[Tuple[Chem.Mol, set]],
):
    for t in library:
        if t.n_reactants == 1:
            if groups.intersection(t.triggers[0]):
                yield predict(t, [mol], groups)
        elif t.n_reactants == 2:
            for pmol, pgroups in partners:
                if groups.intersection(t.triggers[0]) and pgroups.intersection(t.triggers[1]):
                    yield predict(t, [mol, pmol], groups | pgroups)
                if pgroups.intersection(t.triggers[0]) and groups.intersection(t.triggers[1]):
                    yield predict(t, [pmol, mol], groups | pgroups)


def expand(
    smiles: str,
    library: Sequence[Transformation],
    include_minor: bool = False,
    partners: Sequence[MoleculeInput] = (),
) -> List[Tuple[Prediction, str, Tuple[str, ...]]]:
    """Every (prediction, product SMILES, warnings) reachable in one step from ``smiles``.

    ``partners`` are extra reagents (e.g. ``"CC(C)=O"`` for a Grignard addition)
    that let two-component transformations take part.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    groups = {g.group.id for g in detect_groups(mol)}
    partner_info = []
    for p in partners:
        pmol = p if isinstance(p, Chem.Mol) else parse_molecule(p)
        partner_info.append((pmol, {g.group.id for g in detect_groups(pmol)}))
    out = []
    for pred in _candidate_predictions(mol, groups, library, partner_info):
        if pred is None or not pred.feasible:
            continue
        warnings = tuple(c.message for c in pred.cautions if c.severity == "warn")
        products = list(pred.major)
        if pred.exhaustive is not None:
            products = [pred.exhaustive] + [p for p in pred.major if p is not pred.exhaustive]
        if include_minor:
            products += pred.minor
        seen = set()
        for p in products:
            if p.smiles not in seen:
                seen.add(p.smiles)
                out.append((pred, p.smiles, warnings))
    return out


def find_routes(
    start: MoleculeInput,
    target: MoleculeInput,
    max_steps: int = 3,
    max_routes: int = 5,
    include_minor: bool = False,
    library: Optional[Sequence[Transformation]] = None,
    exclude: Iterable[str] = (),
    max_nodes: int = 20000,
    partners: Sequence[MoleculeInput] = (),
) -> List[Route]:
    """Find the shortest sequences of library transformations that turn ``start`` into ``target``.

    One-component transformations (with built-in reagents such as PCC, HBr
    or Ph3P=CH2) are always used. Two-component transformations are used
    only with the co-reactants given in ``partners``. If a step gives several
    fragments (ozonolysis, ester cleavage), any fragment can be carried
    forward. Molecules are compared without stereochemistry.

    Returns at most ``max_routes`` routes of the minimum length found, with
    the routes carrying the fewest chemoselectivity warnings first. Returns an
    empty list if no route exists within ``max_steps``.
    """
    start_smiles = canonical_smiles(parse_molecule(start))
    target_smiles = canonical_smiles(parse_molecule(target))
    target_key = _key(target_smiles)
    excluded = set(exclude)
    lib = [t for t in (library if library is not None else load_transformations()) if t.id not in excluded]
    partner_mols = [parse_molecule(p) for p in partners]
    partner_keys = {_key(canonical_smiles(m)) for m in partner_mols}

    if _key(start_smiles) == target_key:
        return [Route(start_smiles, target_smiles, [])]

    seen = {_key(start_smiles)}
    level: Dict[str, Tuple[str, List[List[RouteStep]]]] = {_key(start_smiles): (start_smiles, [[]])}
    explored = 0

    for _depth in range(max_steps):
        nxt: Dict[str, Tuple[str, List[List[RouteStep]]]] = {}
        for _node_key, (node_smiles, node_paths) in level.items():
            explored += 1
            if explored > max_nodes:
                break
            for pred, product, warnings in expand(node_smiles, lib, include_minor, partner_mols):
                for frag in product.split("."):
                    key = _key(frag)
                    if key in seen or key in partner_keys:
                        continue
                    co = tuple(r for r in pred.reactants if r != node_smiles)
                    step = RouteStep(pred.transformation, node_smiles, frag, product, warnings, co)
                    entry = nxt.setdefault(key, (frag, []))
                    for path in node_paths:
                        if len(entry[1]) >= max_routes:
                            break
                        entry[1].append(path + [step])
        if target_key in nxt:
            paths = sorted(nxt[target_key][1], key=lambda p: sum(len(s.warnings) for s in p))
            return [Route(start_smiles, target_smiles, p) for p in paths[:max_routes]]
        if not nxt or explored > max_nodes:
            break
        seen.update(nxt)
        level = nxt
    return []
