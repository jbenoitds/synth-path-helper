import doctest
import json

import pytest

import synth_path_helper
from synth_path_helper import predictor
from synth_path_helper.cli import main


def run(capsys, *argv):
    code = main(list(argv))
    out, err = capsys.readouterr()
    return code, out, err


def test_analyze_text(capsys):
    code, out, _ = run(capsys, "analyze", "CCO")
    assert code == 0
    assert "Primary alcohol" in out
    assert "PCC" in out and "CC=O" in out


def test_shorthand_without_subcommand(capsys):
    code, out, _ = run(capsys, "cyclohexanol", "--brief")
    assert code == 0 and "O=C1CCCCC1" in out
    assert "Conditions" not in out


def test_analyze_json(capsys):
    code, out, _ = run(capsys, "analyze", "CC(=O)O", "CCO", "--json")
    assert code == 0
    data = json.loads(out)
    assert len(data["components"]) == 2
    assert any(p["transformation"]["id"] == "fischer_esterification" for p in data["predictions"])


def test_analyze_svg(capsys, tmp_path):
    svg = tmp_path / "out.svg"
    code, _, _ = run(capsys, "analyze", "CCO", "--svg", str(svg))
    assert code == 0 and svg.read_text(encoding="utf-8").lstrip().startswith("<?xml")


def test_groups_command(capsys):
    code, out, _ = run(capsys, "groups", "OCC(=O)c1ccccc1")
    assert code == 0 and "Ketone" in out and "Primary alcohol" in out
    code, out, _ = run(capsys, "groups", "CCO", "--json")
    assert json.loads(out)[0]["functional_groups"][0]["id"] == "alcohol_primary"


def test_list_and_show(capsys):
    code, out, _ = run(capsys, "list")
    assert code == 0 and "ox_pcc_primary" in out and "== Oxidation ==" in out
    code, out, _ = run(capsys, "show", "hbr_markovnikov")
    assert code == 0 and json.loads(out)["selectivity"] == "markovnikov"


def test_route_command(capsys, tmp_path):
    svg = tmp_path / "route.svg"
    code, out, _ = run(capsys, "route", "cyclohexanol", "C1CCC2OC2C1", "--svg", str(svg))
    assert code == 0 and "dehydration" in out and "epoxidation" in out
    assert svg.exists()
    code, out, _ = run(capsys, "route", "bromobenzene", "CC(C)(O)c1ccccc1", "--with", "acetone", "--json")
    assert code == 0 and json.loads(out)[0]["steps"][1]["co_reactants"] == ["CC(C)=O"]
    code, out, _ = run(capsys, "route", "CCO", "c1ccccc1", "--max-steps", "1")
    assert code == 3 and "No route found" in out


def test_errors(capsys):
    code, _, err = run(capsys, "analyze", "C1CC")
    assert code == 2 and "Could not parse" in err
    code, _, err = run(capsys, "show", "nope")
    assert code == 2 and "Unknown transformation" in err
    code, _, err = run(capsys, "--library", "missing.json", "analyze", "CCO")
    assert code == 2 and "could not load library" in err
    assert main([]) == 1


def test_custom_library_option(capsys, tmp_path):
    lib = tmp_path / "lib.json"
    lib.write_text(json.dumps({"transformations": [{
        "id": "custom_methylation", "name": "O-methylation", "category": "Substitution",
        "triggers": ["alcohol_primary"], "smarts": "[CX4:1][OX2H:2]>>[C:1][O:2]C",
        "reagents": ["MeI"], "solvents": ["DMF"], "conditions": "rt",
    }]}), encoding="utf-8")
    code, out, _ = run(capsys, "--library", str(lib), "CCO", "--only", "custom_methylation")
    assert code == 0 and "CCOC" in out


@pytest.mark.parametrize("module", [synth_path_helper, predictor])
def test_docstring_examples(module):
    result = doctest.testmod(module)
    assert result.failed == 0
