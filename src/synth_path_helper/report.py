"""Human-readable text reports and SVG depictions."""

from __future__ import annotations

import textwrap
from typing import Iterable, List, Optional, Sequence

from rdkit import Chem
from rdkit.Chem import Draw

from .predictor import Analysis, Prediction, Product
from .routes import Route
from .selectivity import RULE_DESCRIPTIONS
from .transformations import Transformation, categories

WIDTH = 88
LABEL = 15  # width of the "Reagents      :" column

SEVERITY_MARK = {"info": "(i)", "warn": "(!)", "fails": "(x)"}


def _wrap(label: str, text: str, indent: int = 4) -> List[str]:
    prefix = " " * indent + f"{label:<{LABEL}}: "
    body_indent = " " * len(prefix)
    if not text:
        text = "-"
    return textwrap.wrap(text, width=WIDTH, initial_indent=prefix, subsequent_indent=body_indent) or [prefix]


def _join(items: Sequence[str]) -> str:
    return "; ".join(items) if items else "-"


def _product_line(p: Product) -> str:
    return f"{p.smiles}  [{p.formula}, {p.mol_weight:g} g/mol]"


def format_prediction(pred: Prediction, verbose: bool = True) -> List[str]:
    t = pred.transformation
    status = "" if pred.feasible else "  ** NOT EXPECTED TO WORK **"
    lines = [f"  [{t.id}] {t.name}{status}"]
    from .functional_groups import group_name

    if len(pred.reactants) > 1:
        lines += _wrap("Reactants", " + ".join(pred.reactants))
    lines += _wrap("Triggered by", ", ".join(group_name(g) for g in pred.triggered_by))
    lines += _wrap("Reagents", _join(t.reagents))
    lines += _wrap("Catalyst(s)", _join(t.catalysts))
    lines += _wrap("Solvent(s)", _join(t.solvents))
    if verbose:
        lines += _wrap("Conditions", t.conditions)
        if t.mechanism:
            lines += _wrap("Mechanism", t.mechanism)
    label = "Major product" if len(pred.major) == 1 else "Major products"
    for i, p in enumerate(pred.major):
        lines += _wrap(label if i == 0 else "", _product_line(p))
    if pred.exhaustive is not None:
        lines += _wrap("", f"(excess reagent: all {pred.reactive_sites} reactive sites converted)")
    for i, p in enumerate(pred.minor):
        minor_label = "With 1 equiv" if pred.exhaustive is not None else "Minor"
        lines += _wrap(minor_label if i == 0 else "", _product_line(p))
    if t.selectivity and verbose:
        lines += _wrap("Selectivity", RULE_DESCRIPTIONS.get(t.selectivity["rule"], t.selectivity["rule"]))
    if verbose and t.alternatives:
        lines += _wrap("Alternatives", _join(t.alternatives))
    for c in pred.cautions:
        lines += _wrap(SEVERITY_MARK[c.severity] + " Caution", c.message)
    return lines


def format_analysis(analysis: Analysis, verbose: bool = True) -> str:
    out: List[str] = []
    title = f"synth-path-helper :: {analysis.query}"
    out += [title, "=" * min(len(title), WIDTH)]
    for idx, comp in enumerate(analysis.components, 1):
        head = f"Component {idx}: " if len(analysis.components) > 1 else "Input: "
        out.append(f"{head}{comp.smiles}   ({comp.formula}, {comp.mol_weight:g} g/mol)")
        out.append("")
        out.append("  Functional groups")
        if not comp.groups:
            out.append("    (none recognised)")
        for gm in comp.groups:
            out += textwrap.wrap(
                f"{gm.group.name} x{gm.count}: {gm.group.notes}",
                width=WIDTH,
                initial_indent="    - ",
                subsequent_indent="      ",
            )
        out.append("")

    preds = analysis.predictions
    out.append(f"Predicted transformations ({len(preds)})")
    out.append("-" * WIDTH)
    if not preds:
        out.append("  No transformation in the library applies to this input.")
    order = categories()
    by_cat = {}
    for p in preds:
        by_cat.setdefault(p.transformation.category, []).append(p)
    for cat in sorted(by_cat, key=lambda c: order.index(c) if c in order else len(order)):
        out.append(f"\n== {cat} ==")
        for p in by_cat[cat]:
            out += format_prediction(p, verbose)
            out.append("")

    if analysis.partner_hints:
        out.append("Possible with a second reactant")
        out.append("-" * WIDTH)
        for h in analysis.partner_hints:
            out += textwrap.wrap(
                f"[{h.transformation.id}] {h.transformation.name}: {h.message}",
                width=WIDTH, initial_indent="  - ", subsequent_indent="    ",
            )
    return "\n".join(out).rstrip() + "\n"


def format_routes(routes: Sequence[Route], start: str, target: str) -> str:
    out = [f"Route search: {start}  ->  {target}", "=" * WIDTH]
    if not routes:
        out.append("No route found within the step limit. Try --max-steps 4 or --include-minor.")
        return "\n".join(out) + "\n"
    for n, route in enumerate(routes, 1):
        out.append(f"\nRoute {n} ({len(route)} step{'s' if len(route) != 1 else ''})")
        for i, step in enumerate(route.steps, 1):
            t = step.transformation
            lhs = " + ".join((step.reactant,) + step.co_reactants)
            out.append(f"  Step {i}: {lhs}  ->  {step.product}")
            out += _wrap("Reaction", f"{t.name} [{t.id}]", indent=6)
            out += _wrap("Reagents", _join(t.reagents), indent=6)
            if t.catalysts:
                out += _wrap("Catalyst(s)", _join(t.catalysts), indent=6)
            out += _wrap("Solvent(s)", _join(t.solvents), indent=6)
            if step.full_product != step.product:
                out += _wrap("Also formed", step.full_product, indent=6)
            for w in step.warnings:
                out += _wrap("(!) Caution", w, indent=6)
    return "\n".join(out) + "\n"


def format_library(library: Iterable[Transformation]) -> str:
    out = []
    current = None
    for t in sorted(library, key=lambda t: (categories().index(t.category) if t.category in categories() else 99)):
        if t.category != current:
            current = t.category
            out.append(f"\n== {current} ==")
        out.append(f"  {t.id:<26} {t.name}")
    return "\n".join(out).lstrip("\n") + "\n"


def draw_analysis_svg(analysis: Analysis, mols_per_row: int = 4, size=(260, 200)) -> str:
    """An SVG grid: the input followed by the main product of each feasible prediction."""
    mols = [Chem.MolFromSmiles(analysis.smiles)]
    legends = ["input"]
    seen = set()
    for pred in analysis.predictions:
        if not pred.feasible:
            continue
        smi = pred.main_product.smiles
        if smi in seen:
            continue
        seen.add(smi)
        mols.append(Chem.MolFromSmiles(smi))
        legends.append(pred.transformation.id)
    return Draw.MolsToGridImage(
        mols, molsPerRow=mols_per_row, subImgSize=size, legends=legends, useSVG=True
    )


def draw_route_svg(route: Route, size=(260, 200)) -> str:
    mols = [Chem.MolFromSmiles(route.start)] + [Chem.MolFromSmiles(s.product) for s in route.steps]
    legends = ["start"] + [f"{i}. {s.transformation.id}" for i, s in enumerate(route.steps, 1)]
    return Draw.MolsToGridImage(mols, molsPerRow=len(mols), subImgSize=size, legends=legends, useSVG=True)
