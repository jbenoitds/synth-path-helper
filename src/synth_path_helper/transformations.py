"""The transformation library: loading, validating and querying reaction templates."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from rdkit.Chem import AllChem, rdChemReactions

from .chem import DATA_DIR
from .functional_groups import load_groups

SEVERITIES = ("info", "warn", "fails")
SELECTIVITY_RULES = ("markovnikov", "anti_markovnikov", "zaitsev", "hofmann", "eas")


@dataclass(frozen=True)
class Caution:
    """A chemoselectivity note that applies when another group is present."""

    if_present: Tuple[str, ...]
    message: str
    severity: str = "warn"
    always: bool = False  # also apply when the group is one of the template's triggers


@dataclass(frozen=True)
class Transformation:
    id: str
    name: str
    category: str
    smarts: Tuple[str, ...]
    triggers: Tuple[Tuple[str, ...], ...]
    reagents: Tuple[str, ...] = ()
    catalysts: Tuple[str, ...] = ()
    solvents: Tuple[str, ...] = ()
    conditions: str = ""
    mechanism: str = ""
    alternatives: Tuple[str, ...] = ()
    selectivity: Optional[Dict[str, Any]] = field(default=None, compare=False, hash=False)
    exhaustive: bool = False
    cautions: Tuple[Caution, ...] = ()
    partner_hint: str = ""
    example: Optional[Dict[str, str]] = field(default=None, compare=False, hash=False)
    reactions: Tuple[rdChemReactions.ChemicalReaction, ...] = field(
        default=(), compare=False, hash=False, repr=False
    )

    def __post_init__(self) -> None:
        rxns = []
        for sm in self.smarts:
            try:
                rxn = AllChem.ReactionFromSmarts(sm)
            except ValueError as exc:  # RDKit raises a parser exception for malformed input
                raise ValueError(f"Invalid reaction SMARTS in {self.id!r}: {sm} ({exc})") from None
            if rxn is None:
                raise ValueError(f"Invalid reaction SMARTS in {self.id!r}: {sm}")
            rxn.Initialize()
            rxns.append(rxn)
        counts = {r.GetNumReactantTemplates() for r in rxns}
        if len(counts) != 1:
            raise ValueError(f"All templates of {self.id!r} must have the same number of reactants")
        n = counts.pop()
        if n != len(self.triggers):
            raise ValueError(
                f"{self.id!r}: {n} reactant template(s) but {len(self.triggers)} trigger list(s)"
            )
        for rxn in rxns:
            if rxn.GetNumProductTemplates() != 1:
                raise ValueError(
                    f"{self.id!r}: use a single (grouped) product template, e.g. '>>(A.B)'"
                )
        object.__setattr__(self, "reactions", tuple(rxns))

    @property
    def n_reactants(self) -> int:
        return len(self.triggers)

    @property
    def all_triggers(self) -> Tuple[str, ...]:
        return tuple(t for group in self.triggers for t in group)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "category": self.category,
            "mechanism": self.mechanism,
            "reagents": list(self.reagents),
            "catalysts": list(self.catalysts),
            "solvents": list(self.solvents),
            "conditions": self.conditions,
            "alternatives": list(self.alternatives),
            "selectivity": self.selectivity["rule"] if self.selectivity else None,
            "smarts": list(self.smarts),
        }


def _as_tuple(value: Union[None, str, Iterable[str]]) -> Tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    return tuple(value)


def _parse_entry(entry: Dict[str, Any], known_groups: Sequence[str]) -> Transformation:
    tid = entry["id"]
    triggers_raw = entry["triggers"]
    if triggers_raw and isinstance(triggers_raw[0], str):
        triggers: Tuple[Tuple[str, ...], ...] = (tuple(triggers_raw),)
    else:
        triggers = tuple(tuple(t) for t in triggers_raw)
    for gid in (g for grp in triggers for g in grp):
        if gid not in known_groups:
            raise ValueError(f"{tid!r} triggers on unknown functional group {gid!r}")

    cautions = []
    for c in entry.get("cautions", []):
        sev = c.get("severity", "warn")
        if sev not in SEVERITIES:
            raise ValueError(f"{tid!r}: unknown caution severity {sev!r}")
        groups = _as_tuple(c["if_present"])
        for gid in groups:
            if gid not in known_groups:
                raise ValueError(f"{tid!r}: caution refers to unknown group {gid!r}")
        cautions.append(Caution(groups, c["message"], sev, bool(c.get("always", False))))

    selectivity = entry.get("selectivity")
    if selectivity is not None and selectivity.get("rule") not in SELECTIVITY_RULES:
        raise ValueError(f"{tid!r}: unknown selectivity rule {selectivity.get('rule')!r}")

    return Transformation(
        id=tid,
        name=entry["name"],
        category=entry.get("category", "Other"),
        smarts=_as_tuple(entry["smarts"]),
        triggers=triggers,
        reagents=_as_tuple(entry.get("reagents")),
        catalysts=_as_tuple(entry.get("catalysts")),
        solvents=_as_tuple(entry.get("solvents")),
        conditions=entry.get("conditions", ""),
        mechanism=entry.get("mechanism", ""),
        alternatives=_as_tuple(entry.get("alternatives")),
        selectivity=selectivity,
        exhaustive=bool(entry.get("exhaustive", False)),
        cautions=tuple(cautions),
        partner_hint=entry.get("partner_hint", ""),
        example=entry.get("example"),
    )


def load_library(path: Union[str, Path]) -> List[Transformation]:
    """Load and validate a transformation library JSON file."""
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    known = [g.id for g in load_groups()]
    out: List[Transformation] = []
    seen = set()
    for entry in data["transformations"]:
        t = _parse_entry(entry, known)
        if t.id in seen:
            raise ValueError(f"Duplicate transformation id {t.id!r} in {path}")
        seen.add(t.id)
        out.append(t)
    return out


@lru_cache(maxsize=1)
def _builtin() -> Tuple[Transformation, ...]:
    return tuple(load_library(DATA_DIR / "transformations.json"))


def load_transformations(extra_paths: Sequence[Union[str, Path]] = ()) -> List[Transformation]:
    """Return the built-in library plus any user libraries.

    User entries with the same id as a built-in entry replace it.
    """
    library = {t.id: t for t in _builtin()}
    for path in extra_paths:
        for t in load_library(path):
            library[t.id] = t
    return list(library.values())


def get_transformation(tid: str, library: Optional[Sequence[Transformation]] = None) -> Transformation:
    for t in library or _builtin():
        if t.id == tid:
            return t
    raise KeyError(f"Unknown transformation {tid!r}")


def categories(library: Optional[Sequence[Transformation]] = None) -> List[str]:
    seen: List[str] = []
    for t in library or _builtin():
        if t.category not in seen:
            seen.append(t.category)
    return seen
