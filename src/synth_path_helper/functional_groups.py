"""Functional-group perception with SMARTS patterns."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

from rdkit import Chem

from .chem import DATA_DIR


@dataclass(frozen=True)
class FunctionalGroup:
    id: str
    name: str
    family: str
    smarts: str
    notes: str = ""
    pattern: Chem.Mol = field(default=None, compare=False, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        patt = Chem.MolFromSmarts(self.smarts)
        if patt is None:
            raise ValueError(f"Invalid SMARTS for functional group {self.id!r}: {self.smarts}")
        object.__setattr__(self, "pattern", patt)


@dataclass(frozen=True)
class GroupMatch:
    """A functional group found in a molecule."""

    group: FunctionalGroup
    matches: Tuple[Tuple[int, ...], ...]

    @property
    def count(self) -> int:
        return len(self.matches)

    @property
    def anchor_atoms(self) -> List[int]:
        """Index of the first pattern atom in each match (e.g. the OH oxygen)."""
        return [m[0] for m in self.matches]

    def to_dict(self) -> Dict[str, object]:
        return {
            "id": self.group.id,
            "name": self.group.name,
            "family": self.group.family,
            "count": self.count,
            "atoms": [list(m) for m in self.matches],
        }


def _load(path: Union[str, Path]) -> List[FunctionalGroup]:
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    groups = []
    seen: Set[str] = set()
    for entry in data["groups"]:
        if entry["id"] in seen:
            raise ValueError(f"Duplicate functional group id {entry['id']!r} in {path}")
        seen.add(entry["id"])
        groups.append(
            FunctionalGroup(
                id=entry["id"],
                name=entry["name"],
                family=entry.get("family", "Other"),
                smarts=entry["smarts"],
                notes=entry.get("notes", ""),
            )
        )
    return groups


@lru_cache(maxsize=1)
def load_groups() -> Tuple[FunctionalGroup, ...]:
    """Return the built-in functional-group library (cached)."""
    return tuple(_load(DATA_DIR / "functional_groups.json"))


def get_group(group_id: str) -> FunctionalGroup:
    for g in load_groups():
        if g.id == group_id:
            return g
    raise KeyError(f"Unknown functional group {group_id!r}")


def group_name(group_id: str) -> str:
    try:
        return get_group(group_id).name
    except KeyError:
        return group_id


def detect_groups(
    mol: Chem.Mol, groups: Optional[Sequence[FunctionalGroup]] = None
) -> List[GroupMatch]:
    """Find every functional group in ``mol``, in library order."""
    found = []
    for group in groups or load_groups():
        matches = mol.GetSubstructMatches(group.pattern, uniquify=True)
        if matches:
            found.append(GroupMatch(group=group, matches=tuple(matches)))
    return found


def group_ids(mol: Chem.Mol) -> Set[str]:
    return {m.group.id for m in detect_groups(mol)}
