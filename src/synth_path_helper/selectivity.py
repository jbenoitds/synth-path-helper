"""Rule-based ranking of regioisomeric products.

Each rule scores one product molecule; a higher score means a more favoured
product. The product must still carry the ``old_mapno`` atom properties that
RDKit adds in ``RunReactants``, because the rules look at the atoms that were
mapped in the template (``site``, ``other``, ``atoms``).

Rules
-----
markovnikov       The new heteroatom sits on the carbon that would carry the
                  more stable carbocation (more substituted, benzylic).
anti_markovnikov  The reverse (radical HBr, hydroboration).
zaitsev           The more substituted (and conjugated) alkene wins.
hofmann           The less substituted alkene wins (bulky bases).
eas               Directing effects of the substituents already on the
                  benzene ring (ortho/para vs meta, activation strength,
                  steric hindrance at ortho).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem

# Products whose score is within this fraction of the best score (or within
# MIN_TOLERANCE in absolute terms) are reported as major products together.
# That is how ortho/para mixtures and 'no preference' cases are shown.
RELATIVE_TOLERANCE = 0.15
MIN_TOLERANCE = 0.15


def _atom_by_mapno(mol: Chem.Mol, mapno: int) -> Optional[Chem.Atom]:
    for atom in mol.GetAtoms():
        if atom.HasProp("old_mapno") and atom.GetIntProp("old_mapno") == mapno:
            return atom
    return None


def _carbon_neighbours(atom: Chem.Atom, exclude: Sequence[int] = ()) -> int:
    return sum(
        1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 6 and nb.GetIdx() not in exclude
    )


def _is_benzylic_or_allylic(atom: Chem.Atom, exclude: Sequence[int] = ()) -> bool:
    for nb in atom.GetNeighbors():
        if nb.GetIdx() in exclude:
            continue
        if nb.GetIsAromatic():
            return True
        for bond in nb.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtomIdx(nb.GetIdx()) != atom.GetIdx():
                return True
    return False


def _cation_stability(atom: Chem.Atom) -> float:
    """Rough carbocation-stability index: alkyl substitution plus a resonance bonus."""
    score = float(_carbon_neighbours(atom))
    if any(nb.GetIsAromatic() for nb in atom.GetNeighbors()):
        score += 1.5
    return score


def _markovnikov(mol: Chem.Mol, spec: Dict[str, Any]) -> float:
    site = _atom_by_mapno(mol, int(spec.get("site", 2)))
    other = _atom_by_mapno(mol, int(spec.get("other", 1)))
    if site is None or other is None:
        return 0.0
    return _cation_stability(site) - _cation_stability(other)


def _alkene_substitution(mol: Chem.Mol, spec: Dict[str, Any]) -> Tuple[float, float]:
    a_no, b_no = spec.get("atoms", [1, 2])
    a, b = _atom_by_mapno(mol, int(a_no)), _atom_by_mapno(mol, int(b_no))
    if a is None or b is None:
        return 0.0, 0.0
    subst = _carbon_neighbours(a, [b.GetIdx()]) + _carbon_neighbours(b, [a.GetIdx()])
    conj = 0.5 if (
        _is_benzylic_or_allylic(a, [b.GetIdx()]) or _is_benzylic_or_allylic(b, [a.GetIdx()])
    ) else 0.0
    return float(subst), conj


def _zaitsev(mol: Chem.Mol, spec: Dict[str, Any]) -> float:
    subst, conj = _alkene_substitution(mol, spec)
    return subst + conj


def _hofmann(mol: Chem.Mol, spec: Dict[str, Any]) -> float:
    subst, _ = _alkene_substitution(mol, spec)
    return -subst


# --- electrophilic aromatic substitution ------------------------------------

ORTHO, META, PARA = 1, 2, 3


def classify_substituent(mol: Chem.Mol, ring_atom: Chem.Atom, sub: Chem.Atom) -> Tuple[str, float, bool]:
    """Classify the substituent attached to ``ring_atom`` through atom ``sub``.

    Returns ``(director, weight, bulky)``: director is ``"op"`` (ortho/para)
    or ``"meta"``; weight is how strongly it directs; ``bulky`` makes the
    ortho positions less favourable.
    """
    z = sub.GetAtomicNum()
    charge = sub.GetFormalCharge()

    def double_bonded_heteroatom(atom: Chem.Atom) -> bool:
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetIdx() == ring_atom.GetIdx():
                continue
            if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) and other.GetAtomicNum() in (7, 8, 16):
                return True
        return False

    if z in (9, 17, 35, 53):
        return "op", 1.0, False  # deactivating, but an o/p director
    if z == 8:
        if charge < 0:
            return "op", 3.0, False
        acyl = any(
            nb.GetIdx() != ring_atom.GetIdx() and double_bonded_heteroatom(nb) for nb in sub.GetNeighbors()
        )
        return ("op", 2.0, False) if acyl else ("op", 3.0, False)
    if z == 7:
        if charge > 0:
            return "meta", 3.0, False  # -NO2, -NR3+
        acyl = any(
            nb.GetIdx() != ring_atom.GetIdx() and double_bonded_heteroatom(nb) for nb in sub.GetNeighbors()
        )
        bulky = sub.GetDegree() >= 3
        return ("op", 2.0, bulky) if acyl else ("op", 3.0, bulky)
    if z == 16:
        if double_bonded_heteroatom(sub):
            return "meta", 2.0, False  # -SO3H, -SO2R
        return "op", 1.0, False
    if z == 6:
        if sub.GetIsAromatic():
            return "op", 1.0, False  # phenyl or ring fusion
        if double_bonded_heteroatom(sub):
            return "meta", 2.0, False  # C=O, C#N
        halogens = sum(1 for nb in sub.GetNeighbors() if nb.GetAtomicNum() in (9, 17))
        if halogens >= 2:
            return "meta", 2.0, False  # CF3, CCl3
        bulky = sub.GetDegree() >= 4  # tert-butyl-like quaternary carbon
        return "op", 1.0, bulky
    return "op", 0.5, False


def _position_score(director: str, position: int, bulky: bool) -> float:
    if director == "meta":
        return 2.0 if position == META else 0.0
    if position == PARA:
        return 2.0
    if position == ORTHO:
        return 0.6 if bulky else 1.75
    return 0.0


def eas_breakdown(mol: Chem.Mol, site_idx: int) -> List[Dict[str, Any]]:
    """List the substituents on the ring that contains ``site_idx`` and how each one directs."""
    ring = None
    for r in mol.GetRingInfo().AtomRings():
        if site_idx in r and len(r) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
            ring = r
            break
    if ring is None:
        return []
    dmat = Chem.GetDistanceMatrix(mol)
    out = []
    for ring_idx in ring:
        if ring_idx == site_idx:
            continue
        ring_atom = mol.GetAtomWithIdx(ring_idx)
        for nb in ring_atom.GetNeighbors():
            if nb.GetIdx() in ring:
                continue
            director, weight, bulky = classify_substituent(mol, ring_atom, nb)
            pos = int(dmat[ring_idx][site_idx])
            out.append(
                {
                    "ring_atom": ring_idx,
                    "substituent_atom": nb.GetIdx(),
                    "director": director,
                    "weight": weight,
                    "bulky": bulky,
                    "position": {ORTHO: "ortho", META: "meta", PARA: "para"}.get(pos, "?"),
                    # ring fusion (naphthalene) rather than a separate substituent (biphenyl)
                    "fused": mol.GetBondBetweenAtoms(ring_idx, nb.GetIdx()).IsInRing(),
                    "_pos": pos,
                }
            )
    return out


def _eas(mol: Chem.Mol, spec: Dict[str, Any]) -> float:
    site = _atom_by_mapno(mol, int(spec.get("site", 1)))
    if site is None:
        return 0.0
    subs = eas_breakdown(mol, site.GetIdx())
    if not subs:
        return 0.0  # unsubstituted benzene: every position is equivalent
    score = 0.0
    ortho_count = 0
    for s in subs:
        score += s["weight"] * _position_score(s["director"], s["_pos"], s["bulky"])
        if s["_pos"] == ORTHO:
            ortho_count += 1
        # Alpha position of a fused ring system (naphthalene C1) is kinetically favoured.
        if s["fused"] and s["_pos"] == ORTHO:
            score += 1.0
    if ortho_count >= 2:
        score -= 1.5  # squeezed between two groups
    return score


_RULES = {
    "markovnikov": _markovnikov,
    "anti_markovnikov": lambda m, s: -_markovnikov(m, s),
    "zaitsev": _zaitsev,
    "hofmann": _hofmann,
    "eas": _eas,
}

RULE_DESCRIPTIONS = {
    "markovnikov": "Markovnikov: the new group ends up on the carbon that best stabilises a positive charge.",
    "anti_markovnikov": "Anti-Markovnikov: the new group ends up on the less substituted carbon.",
    "zaitsev": "Zaitsev: the more substituted (more stable) alkene is the major product.",
    "hofmann": "Hofmann: a bulky base gives the less substituted alkene.",
    "eas": "EAS directing effects: activators and halogens direct ortho/para, deactivators direct meta.",
}


def score_product(mol: Chem.Mol, spec: Optional[Dict[str, Any]]) -> float:
    """Score ``mol`` with the selectivity rule ``spec`` (0.0 when there is no rule)."""
    if not spec:
        return 0.0
    return float(_RULES[spec["rule"]](mol, spec))


def is_major(score: float, best: float) -> bool:
    tolerance = max(RELATIVE_TOLERANCE * abs(best), MIN_TOLERANCE)
    return score >= best - tolerance
