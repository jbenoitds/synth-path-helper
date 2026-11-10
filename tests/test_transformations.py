import json

import pytest

from synth_path_helper import analyze, get_transformation, load_transformations
from synth_path_helper.chem import canonicalize
from synth_path_helper.transformations import SEVERITIES, categories, load_library

LIBRARY = load_transformations()


def test_library_size_and_unique_ids():
    assert len(LIBRARY) >= 60
    assert len({t.id for t in LIBRARY}) == len(LIBRARY)


@pytest.mark.parametrize("t", LIBRARY, ids=lambda t: t.id)
def test_entry_is_complete(t):
    assert t.name and t.category and t.conditions
    assert t.reagents, "every transformation must name its reagents"
    assert t.solvents, "every transformation must name a solvent"
    assert t.example and "reactants" in t.example and "major" in t.example
    for c in t.cautions:
        assert c.severity in SEVERITIES
        assert "{group}" in c.message
    if t.n_reactants == 2:
        assert t.partner_hint, "two-component reactions need a partner hint"


@pytest.mark.parametrize("t", LIBRARY, ids=lambda t: t.id)
def test_documented_example(t):
    """Every library entry predicts its own documented example."""
    report = analyze(t.example["reactants"], ids=[t.id])
    preds = report.find(t.id)
    assert preds, f"{t.id} did not fire on {t.example['reactants']}"
    expected = canonicalize(t.example["major"])
    majors = [p.smiles for pred in preds for p in pred.major]
    assert expected in majors


def test_categories_are_ordered_and_known():
    cats = categories()
    assert "Oxidation" in cats and "Reduction" in cats and "Aromatic substitution" in cats


def test_get_transformation():
    t = get_transformation("ox_pcc_primary")
    assert "PCC" in t.reagents[0]
    assert t.to_dict()["id"] == "ox_pcc_primary"
    with pytest.raises(KeyError):
        get_transformation("nope")


def _write_library(tmp_path, entries):
    path = tmp_path / "lib.json"
    path.write_text(json.dumps({"transformations": entries}), encoding="utf-8")
    return path


BASE_ENTRY = {
    "id": "custom_methylation",
    "name": "O-methylation with MeI / Ag2O",
    "category": "Substitution",
    "triggers": ["alcohol_primary"],
    "smarts": "[CX4:1][OX2H:2]>>[C:1][O:2]C",
    "reagents": ["MeI", "Ag2O"],
    "solvents": ["DMF"],
    "conditions": "rt",
    "example": {"reactants": "CCO", "major": "CCOC"},
}


def test_user_library_extends_builtin(tmp_path):
    path = _write_library(tmp_path, [BASE_ENTRY])
    lib = load_transformations([path])
    assert len(lib) == len(LIBRARY) + 1
    report = analyze("CCO", library=lib, ids=["custom_methylation"])
    assert report.find("custom_methylation")[0].major[0].smiles == "CCOC"


def test_user_library_overrides_builtin(tmp_path):
    entry = dict(BASE_ENTRY, id="ox_pcc_primary")
    lib = load_transformations([_write_library(tmp_path, [entry])])
    assert len(lib) == len(LIBRARY)
    assert get_transformation("ox_pcc_primary", lib).reagents == ("MeI", "Ag2O")


@pytest.mark.parametrize(
    "patch, message",
    [
        ({"smarts": "this is not smarts>>"}, "Invalid reaction SMARTS"),
        ({"triggers": ["unicorn"]}, "unknown functional group"),
        ({"triggers": [["alcohol_primary"], ["ketone"]]}, "trigger list"),
        ({"smarts": "[C:1][OX2H:2]>>[C:1]=O.[O:2]"}, r"single \(grouped\) product"),
        ({"selectivity": {"rule": "magic"}}, "unknown selectivity rule"),
        ({"cautions": [{"if_present": "ketone", "message": "x", "severity": "panic"}]}, "severity"),
    ],
)
def test_invalid_entries_are_rejected(tmp_path, patch, message):
    path = _write_library(tmp_path, [dict(BASE_ENTRY, **patch)])
    with pytest.raises(ValueError, match=message):
        load_library(path)


def test_duplicate_ids_rejected(tmp_path):
    path = _write_library(tmp_path, [BASE_ENTRY, BASE_ENTRY])
    with pytest.raises(ValueError, match="Duplicate"):
        load_library(path)
