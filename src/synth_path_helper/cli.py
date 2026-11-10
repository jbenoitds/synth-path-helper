"""Command-line interface.

Examples
--------
    synth-path-helper analyze "CCO"
    synth-path-helper analyze cyclohexanol --category Oxidation
    synth-path-helper analyze "CC(=O)O.CCO" --json
    synth-path-helper route 1-propanol "CCC(=O)O"
    synth-path-helper groups "OCC(=O)c1ccccc1"
    synth-path-helper list
    synth-path-helper show ox_pcc_primary
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import List, Optional

from . import __version__
from .chem import InvalidMoleculeError, canonical_smiles, parse_mixture, parse_molecule
from .functional_groups import detect_groups
from .predictor import analyze
from .report import (
    draw_analysis_svg,
    draw_route_svg,
    format_analysis,
    format_library,
    format_routes,
)
from .routes import find_routes
from .transformations import get_transformation, load_transformations

COMMANDS = ("analyze", "route", "groups", "list", "show")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="synth-path-helper",
        description="Rule-based functional-group analysis and reaction prediction from SMILES.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--library", action="append", default=[], metavar="JSON",
        help="extra transformation library to load (repeatable; same ids override built-ins)",
    )
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("analyze", help="detect functional groups and predict products")
    p.add_argument("molecules", nargs="+", help="SMILES (dot-separated for mixtures) or a common name")
    p.add_argument("--category", action="append", help="only show this category (repeatable)")
    p.add_argument("--only", action="append", metavar="ID", help="only run this transformation id (repeatable)")
    p.add_argument("--feasible", action="store_true", help="hide transformations flagged as not expected to work")
    p.add_argument("--brief", action="store_true", help="shorter output (no conditions or alternatives)")
    p.add_argument("--json", action="store_true", help="print machine-readable JSON")
    p.add_argument("--svg", metavar="FILE", help="also write an SVG grid of the input and products")

    r = sub.add_parser("route", help="search for a multi-step route from START to TARGET")
    r.add_argument("start")
    r.add_argument("target")
    r.add_argument("--max-steps", type=int, default=3)
    r.add_argument("--max-routes", type=int, default=5)
    r.add_argument("--include-minor", action="store_true", help="also follow minor products")
    r.add_argument("--exclude", action="append", default=[], metavar="ID", help="do not use this transformation")
    r.add_argument(
        "--with", dest="partners", action="append", default=[], metavar="SMILES",
        help="co-reactant available for two-component steps, e.g. --with 'CC(C)=O' (repeatable)",
    )
    r.add_argument("--json", action="store_true")
    r.add_argument("--svg", metavar="FILE", help="write an SVG of the first route")

    g = sub.add_parser("groups", help="list the functional groups in a molecule")
    g.add_argument("molecules", nargs="+")
    g.add_argument("--json", action="store_true")

    sub.add_parser("list", help="list every transformation in the library")

    s = sub.add_parser("show", help="show the details of one transformation")
    s.add_argument("id")
    return parser


def _write(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # 'synth-path-helper CCO' is shorthand for 'synth-path-helper analyze CCO'
    idx = 0
    while idx < len(argv):
        arg = argv[idx]
        if arg == "--library":
            idx += 2  # skip the option's value
            continue
        if arg.startswith("-"):
            idx += 1
            continue
        if arg not in COMMANDS:
            # keep global options such as --library before the inserted command
            argv.insert(idx, "analyze")
        break
    # Windows consoles may not encode Greek letters such as the 'ν' in 'hν'
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")  # type: ignore[attr-defined]
        except (AttributeError, ValueError):
            pass

    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 1

    try:
        library = load_transformations(args.library)
    except (OSError, ValueError, KeyError) as exc:
        print(f"error: could not load library: {exc}", file=sys.stderr)
        return 2

    try:
        if args.command == "analyze":
            result = analyze(
                args.molecules if len(args.molecules) > 1 else args.molecules[0],
                library=library,
                categories=args.category,
                ids=args.only,
                feasible_only=args.feasible,
            )
            if args.json:
                print(json.dumps(result.to_dict(), indent=2))
            else:
                print(format_analysis(result, verbose=not args.brief), end="")
            if args.svg:
                _write(args.svg, draw_analysis_svg(result))
                print(f"SVG written to {args.svg}", file=sys.stderr)

        elif args.command == "route":
            routes = find_routes(
                args.start, args.target,
                max_steps=args.max_steps, max_routes=args.max_routes,
                include_minor=args.include_minor, library=library, exclude=args.exclude,
                partners=args.partners,
            )
            start = canonical_smiles(parse_molecule(args.start))
            target = canonical_smiles(parse_molecule(args.target))
            if args.json:
                print(json.dumps([r.to_dict() for r in routes], indent=2))
            else:
                print(format_routes(routes, start, target), end="")
            if args.svg and routes:
                _write(args.svg, draw_route_svg(routes[0]))
                print(f"SVG written to {args.svg}", file=sys.stderr)
            return 0 if routes else 3

        elif args.command == "groups":
            comps = parse_mixture(args.molecules)
            data = [
                {"smiles": canonical_smiles(m), "functional_groups": [g.to_dict() for g in detect_groups(m)]}
                for m in comps
            ]
            if args.json:
                print(json.dumps(data, indent=2))
            else:
                for entry in data:
                    print(entry["smiles"])
                    if not entry["functional_groups"]:
                        print("  (none recognised)")
                    for g in entry["functional_groups"]:
                        atoms = ", ".join(str(a[0]) for a in g["atoms"])
                        print(f"  - {g['name']} x{g['count']}  (anchor atom idx: {atoms})")

        elif args.command == "list":
            print(format_library(library), end="")

        elif args.command == "show":
            t = get_transformation(args.id, library)
            print(json.dumps({**t.to_dict(), "example": t.example,
                              "cautions": [c.__dict__ for c in t.cautions]}, indent=2, default=list))
    except InvalidMoleculeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyError as exc:
        print(f"error: {exc.args[0] if exc.args else exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
