"""Match molecules against the transformation library and predict products."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

from rdkit import Chem

from . import selectivity
from .chem import (
    MoleculeInput,
    canonical_smiles,
    finalize_product,
    formula,
    mol_weight,
    parse_mixture,
)
from .functional_groups import GroupMatch, detect_groups, group_name
from .transformations import Transformation, load_transformations

MAX_EXHAUSTIVE_ITERATIONS = 25
MAX_PRODUCTS_PER_TEMPLATE = 500


@dataclass(frozen=True)
class Product:
    smiles: str
    formula: str
    mol_weight: float
    score: float = 0.0

    @property
    def fragments(self) -> List[str]:
        return self.smiles.split(".")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "formula": self.formula,
            "mol_weight": self.mol_weight,
            "score": round(self.score, 3),
        }


@dataclass(frozen=True)
class CautionHit:
    severity: str
    message: str
    group: str

    def to_dict(self) -> Dict[str, str]:
        return {"severity": self.severity, "message": self.message, "group": self.group}


@dataclass
class Prediction:
    """The predicted outcome of one transformation on one set of reactants."""

    transformation: Transformation
    reactants: List[str]
    triggered_by: List[str]
    major: List[Product]
    minor: List[Product] = field(default_factory=list)
    exhaustive: Optional[Product] = None
    reactive_sites: int = 1
    cautions: List[CautionHit] = field(default_factory=list)

    @property
    def feasible(self) -> bool:
        return not any(c.severity == "fails" for c in self.cautions)

    @property
    def main_product(self) -> Product:
        """The single product to carry forward: the exhaustive product if all sites react, else the top major product."""
        return self.exhaustive or self.major[0]

    def to_dict(self) -> Dict[str, Any]:
        t = self.transformation
        return {
            "transformation": t.to_dict(),
            "reactants": self.reactants,
            "triggered_by": [{"id": g, "name": group_name(g)} for g in self.triggered_by],
            "major_products": [p.to_dict() for p in self.major],
            "minor_products": [p.to_dict() for p in self.minor],
            "exhaustive_product": self.exhaustive.to_dict() if self.exhaustive else None,
            "reactive_sites": self.reactive_sites,
            "feasible": self.feasible,
            "cautions": [c.to_dict() for c in self.cautions],
            "selectivity_rule": selectivity.RULE_DESCRIPTIONS.get(t.selectivity["rule"]) if t.selectivity else None,
        }


@dataclass
class Component:
    smiles: str
    formula: str
    mol_weight: float
    groups: List[GroupMatch]

    @property
    def group_ids(self) -> Set[str]:
        return {g.group.id for g in self.groups}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "formula": self.formula,
            "mol_weight": self.mol_weight,
            "functional_groups": [g.to_dict() for g in self.groups],
        }


@dataclass
class PartnerHint:
    transformation: Transformation
    component: str
    message: str

    def to_dict(self) -> Dict[str, str]:
        return {"transformation": self.transformation.id, "name": self.transformation.name,
                "component": self.component, "message": self.message}


@dataclass
class Analysis:
    query: str
    components: List[Component]
    predictions: List[Prediction]
    partner_hints: List[PartnerHint] = field(default_factory=list)

    @property
    def smiles(self) -> str:
        return ".".join(c.smiles for c in self.components)

    def for_group(self, group_id: str) -> List[Prediction]:
        return [p for p in self.predictions if group_id in p.triggered_by]

    def find(self, transformation_id: str) -> List[Prediction]:
        return [p for p in self.predictions if p.transformation.id == transformation_id]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "query": self.query,
            "smiles": self.smiles,
            "components": [c.to_dict() for c in self.components],
            "predictions": [p.to_dict() for p in self.predictions],
            "partner_hints": [h.to_dict() for h in self.partner_hints],
        }


# --------------------------------------------------------------------------
# Core template application
# --------------------------------------------------------------------------


def _make_product(mol: Chem.Mol, score: float = 0.0) -> Product:
    return Product(canonical_smiles(mol), formula(mol), mol_weight(mol), score)


def run_template(t: Transformation, reactants: Sequence[Chem.Mol]) -> Dict[str, Tuple[Chem.Mol, float]]:
    """Apply every template of ``t`` to ``reactants``.

    Returns ``{canonical_smiles: (product_mol, score)}`` with one entry per
    distinct product, keeping the best selectivity score for each.
    """
    results: Dict[str, Tuple[Chem.Mol, float]] = {}
    for rxn in t.reactions:
        outcomes = rxn.RunReactants(tuple(reactants), MAX_PRODUCTS_PER_TEMPLATE)
        for outcome in outcomes:
            raw = outcome[0]
            prod = finalize_product(raw)
            if prod is None:
                continue
            smi = canonical_smiles(prod)
            score = selectivity.score_product(prod, t.selectivity)
            if smi not in results or score > results[smi][1]:
                results[smi] = (prod, score)
    return results


def count_sites(t: Transformation, mol: Chem.Mol) -> int:
    """How many distinct places in ``mol`` match the (first) reactant template."""
    sites = set()
    for rxn in t.reactions:
        patt = rxn.GetReactantTemplate(0)
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            sites.add(frozenset(match))
    return len(sites)


def _exhaustive(t: Transformation, mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Apply a one-component template repeatedly until no site is left (excess reagent)."""
    current = mol
    changed = False
    for _ in range(MAX_EXHAUSTIVE_ITERATIONS):
        nxt = None
        for rxn in t.reactions:
            for outcome in rxn.RunReactants((current,), 10):
                prod = finalize_product(outcome[0])
                if prod is not None:
                    nxt = Chem.MolFromSmiles(canonical_smiles(prod))
                    break
            if nxt is not None:
                break
        if nxt is None:
            break
        current, changed = nxt, True
    return current if changed else None


def _cautions(t: Transformation, present: Set[str]) -> List[CautionHit]:
    triggers = set(t.all_triggers)
    hits = []
    for c in t.cautions:
        for gid in c.if_present:
            if gid in present and (c.always or gid not in triggers):
                hits.append(CautionHit(c.severity, c.message.format(group=group_name(gid)), gid))
    return hits


def _rank(t: Transformation, results: Dict[str, Tuple[Chem.Mol, float]]) -> Tuple[List[Product], List[Product]]:
    products = sorted(
        (_make_product(Chem.MolFromSmiles(smi), score) for smi, (_, score) in results.items()),
        key=lambda p: (-p.score, p.smiles),
    )
    if not t.selectivity:
        return products, []
    best = products[0].score
    major = [p for p in products if selectivity.is_major(p.score, best)]
    minor = [p for p in products if not selectivity.is_major(p.score, best)]
    return major, minor


def predict(
    t: Transformation,
    reactants: Sequence[MoleculeInput],
    present_groups: Optional[Set[str]] = None,
) -> Optional[Prediction]:
    """Predict the outcome of transformation ``t`` on ``reactants``.

    ``reactants`` must contain one molecule per reactant template, in template
    order. Returns ``None`` when the template does not apply.
    """
    mols = [m if isinstance(m, Chem.Mol) else parse_mixture(m)[0] for m in reactants]
    if len(mols) != t.n_reactants:
        raise ValueError(f"{t.id} needs {t.n_reactants} reactant(s), got {len(mols)}")
    if present_groups is None:
        present_groups = set()
        for m in mols:
            present_groups |= {g.group.id for g in detect_groups(m)}

    results = run_template(t, mols)
    reactant_smiles = [canonical_smiles(m) for m in mols]
    # Discard 'products' identical to the starting material.
    for smi in list(results):
        if smi == ".".join(reactant_smiles):
            del results[smi]
    if not results:
        return None

    major, minor = _rank(t, results)
    sites = count_sites(t, mols[0]) if t.n_reactants == 1 else 1

    exhaustive_product = None
    if t.exhaustive and t.n_reactants == 1 and sites > 1 and not t.selectivity:
        full = _exhaustive(t, mols[0])
        if full is not None:
            full_smi = canonical_smiles(full)
            if full_smi not in {p.smiles for p in major}:
                exhaustive_product = _make_product(full)
                # With excess reagent every site reacts; the single-site
                # products are what one equivalent gives.
                minor = major + minor
                major = [exhaustive_product]

    triggered = [g for g in t.all_triggers if g in present_groups]
    return Prediction(
        transformation=t,
        reactants=reactant_smiles,
        triggered_by=list(dict.fromkeys(triggered)),
        major=major,
        minor=minor,
        exhaustive=exhaustive_product,
        reactive_sites=sites,
        cautions=_cautions(t, present_groups),
    )


# --------------------------------------------------------------------------
# High-level API
# --------------------------------------------------------------------------


def _filter_library(
    library: Optional[Sequence[Transformation]],
    categories: Optional[Iterable[str]],
    ids: Optional[Iterable[str]],
) -> List[Transformation]:
    lib = list(library) if library is not None else load_transformations()
    if categories:
        wanted = {c.lower() for c in categories}
        lib = [t for t in lib if t.category.lower() in wanted]
    if ids:
        wanted_ids = set(ids)
        lib = [t for t in lib if t.id in wanted_ids]
    return lib


def analyze(
    query: Union[MoleculeInput, Sequence[MoleculeInput]],
    library: Optional[Sequence[Transformation]] = None,
    categories: Optional[Iterable[str]] = None,
    ids: Optional[Iterable[str]] = None,
    feasible_only: bool = False,
) -> Analysis:
    """Analyse one molecule (or a mixture) and predict every applicable transformation.

    >>> report = analyze("CCO")
    >>> [p.transformation.id for p in report.predictions][:2]
    ['ox_pcc_primary', 'ox_jones_primary']
    >>> report.find("ox_pcc_primary")[0].major[0].smiles
    'CC=O'
    """
    mols = parse_mixture(query)
    lib = _filter_library(library, categories, ids)

    components = [
        Component(canonical_smiles(m), formula(m), mol_weight(m), detect_groups(m)) for m in mols
    ]
    comp_groups = [c.group_ids for c in components]

    predictions: List[Prediction] = []
    hints: List[PartnerHint] = []

    for t in lib:
        if t.n_reactants == 1:
            for mol, groups in zip(mols, comp_groups):
                if not groups.intersection(t.triggers[0]):
                    continue
                pred = predict(t, [mol], groups)
                if pred is not None:
                    predictions.append(pred)
        elif t.n_reactants == 2:
            applied = False
            for i, j in permutations(range(len(mols)), 2):
                gi, gj = comp_groups[i], comp_groups[j]
                if gi.intersection(t.triggers[0]) and gj.intersection(t.triggers[1]):
                    pred = predict(t, [mols[i], mols[j]], gi | gj)
                    if pred is not None:
                        predictions.append(pred)
                        applied = True
            if not applied and t.partner_hint:
                for comp, groups in zip(components, comp_groups):
                    if groups.intersection(t.triggers[0]) or groups.intersection(t.triggers[1]):
                        hints.append(PartnerHint(t, comp.smiles, t.partner_hint))
                        break

    if feasible_only:
        predictions = [p for p in predictions if p.feasible]
    query_text = query if isinstance(query, str) else ".".join(c.smiles for c in components)
    return Analysis(query=query_text, components=components, predictions=predictions, partner_hints=hints)


def predict_products(smiles: MoleculeInput, transformation_id: str) -> List[str]:
    """Shortcut: the major product SMILES of one named transformation."""
    from .transformations import get_transformation

    t = get_transformation(transformation_id)
    mols = parse_mixture(smiles)
    if len(mols) != t.n_reactants:
        raise ValueError(f"{t.id} needs {t.n_reactants} component(s), got {len(mols)}")
    # Components may be given in any order; try each assignment to the templates.
    for ordering in permutations(mols):
        pred = predict(t, list(ordering))
        if pred is not None:
            return [p.smiles for p in pred.major]
    return []
