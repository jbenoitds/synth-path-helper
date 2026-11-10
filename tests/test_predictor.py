import json

import pytest
from rdkit import Chem

from synth_path_helper import InvalidMoleculeError, analyze, parse_mixture, parse_molecule, predict_products
from synth_path_helper.chem import canonicalize


def ids(report):
    return {p.transformation.id for p in report.predictions}


def test_primary_alcohol_oxidation_pathways():
    report = analyze("CCO")
    assert report.find("ox_pcc_primary")[0].major[0].smiles == "CC=O"
    assert report.find("ox_jones_primary")[0].major[0].smiles == "CC(=O)O"
    pcc = report.find("ox_pcc_primary")[0]
    assert "PCC" in pcc.transformation.reagents[0]
    assert pcc.triggered_by == ["alcohol_primary"]
    assert "ox_secondary" not in ids(report)


def test_secondary_alcohol_gives_ketone_not_acid():
    report = analyze("cyclohexanol")
    assert report.find("ox_secondary")[0].major[0].smiles == "O=C1CCCCC1"
    assert "ox_jones_primary" not in ids(report)


def test_tertiary_alcohol_is_not_oxidised():
    report = analyze("CC(C)(C)O")
    assert not [p for p in report.predictions if p.transformation.category == "Oxidation"]
    assert "hx_tertiary" in ids(report)


def test_no_reaction_for_alkane():
    assert analyze("CC").predictions == []


def test_exhaustive_reaction_on_diol():
    pred = analyze("OCCCCO", ids=["ox_pcc_primary"]).predictions[0]
    assert pred.reactive_sites == 2
    assert pred.exhaustive is not None
    assert pred.major[0].smiles == "O=CCCC=O"
    assert "O=CCCCO" in {p.smiles for p in pred.minor}
    assert pred.main_product.smiles == "O=CCCC=O"


def test_diene_hydrogenation_goes_to_completion():
    pred = analyze("C=CC=C", ids=["hydrogenation_alkene"]).predictions[0]
    assert pred.major[0].smiles == "CCCC"


def test_chemoselectivity_cautions():
    keto_ester = "CCOC(=O)CCC(C)=O"
    report = analyze(keto_ester)
    nabh4 = report.find("red_nabh4_ketone")[0]
    assert nabh4.major[0].smiles == canonicalize("CCOC(=O)CCC(C)O")
    assert any(c.severity == "info" and c.group == "ester" for c in nabh4.cautions)
    lialh4 = report.find("red_lialh4_ester")[0]
    assert any(c.severity == "warn" and c.group == "ketone" for c in lialh4.cautions)
    # the ester's own group is not flagged as a side reaction
    assert not any(c.group == "ester" for c in lialh4.cautions)


def test_infeasible_reactions_are_flagged():
    grignard = analyze("OCCCBr").find("grignard_formation")[0]
    assert not grignard.feasible
    assert any(c.severity == "fails" for c in grignard.cautions)

    fc = analyze("O=[N+]([O-])c1ccccc1").find("fc_acylation")[0]
    assert not fc.feasible

    report = analyze("O=[N+]([O-])c1ccccc1", feasible_only=True)
    assert "fc_acylation" not in ids(report)
    assert "eas_nitration" in ids(report)


def test_trigger_group_caution_with_always_flag():
    pred = analyze("CC(Br)CC").find("sn2_hydroxide")[0]
    assert any("E2" in c.message for c in pred.cautions)


def test_two_component_reaction_in_mixture():
    report = analyze("CC(=O)O.CCO")
    fischer = report.find("fischer_esterification")
    assert fischer and fischer[0].major[0].smiles == "CCOC(C)=O"
    assert fischer[0].reactants == ["CC(=O)O", "CCO"]
    # one-component chemistry is still predicted for each component
    assert "ox_pcc_primary" in ids(report) and "acid_to_acyl_chloride" in ids(report)
    assert not report.partner_hints or all(h.transformation.id != "fischer_esterification"
                                           for h in report.partner_hints)


def test_component_order_does_not_matter():
    a = analyze("CCO.CC(=O)O").find("fischer_esterification")[0].major[0].smiles
    b = analyze(["CC(=O)O", "CCO"]).find("fischer_esterification")[0].major[0].smiles
    assert a == b == "CCOC(C)=O"
    assert predict_products("CC(C)=O.CC[Mg]Br", "grignard_carbonyl") == ["CCC(C)(C)O"]


def test_partner_hint_when_second_reactant_missing():
    report = analyze("CC(=O)O")
    hinted = {h.transformation.id for h in report.partner_hints}
    assert "fischer_esterification" in hinted


def test_name_lookup_and_explicit_hydrogens():
    assert analyze("Cyclohexanol").smiles == "OC1CCCCC1"
    # explicit [CH2]/[OH] must not leave stale hydrogen counts on products
    assert predict_products("[CH3][CH2][OH]", "ox_pcc_primary") == ["CC=O"]


def test_nitro_reduction_neutralises_nitrogen():
    assert predict_products("O=[N+]([O-])c1ccccc1", "red_nitro") == ["Nc1ccccc1"]


def test_stereospecific_alkyne_reductions():
    cis = Chem.MolFromSmiles(predict_products("CC#CCC", "lindlar")[0])
    trans = Chem.MolFromSmiles(predict_products("CC#CCC", "dissolving_metal_alkyne")[0])
    def stereo(m):
        return [b.GetStereo() for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE][0]
    assert stereo(cis) == Chem.BondStereo.STEREOZ
    assert stereo(trans) == Chem.BondStereo.STEREOE
    assert predict_products("C#CCC", "lindlar") == ["C=CCC"]


def test_reaction_centre_chirality_removed_but_other_centres_kept():
    # (S)-2-butanol -> 2-butanone: the stereocentre is destroyed
    assert predict_products("C[C@H](O)CC", "ox_secondary") == ["CCC(C)=O"]
    # remote stereocentre is untouched by the oxidation at C1
    product = predict_products("OC[C@@H](C)CC", "ox_pcc_primary")[0]
    assert "@" in product


def test_every_product_is_valid_and_json_serialisable():
    for smiles in ["OCC(=O)c1ccccc1", "C=CC(=O)OCC", "Nc1ccc(Br)cc1", "CC1(C)CO1", "C#Cc1ccccc1"]:
        report = analyze(smiles)
        assert report.predictions
        for pred in report.predictions:
            for p in pred.major + pred.minor:
                assert Chem.MolFromSmiles(p.smiles) is not None, (pred.transformation.id, p.smiles)
        json.dumps(report.to_dict())


def test_category_and_id_filters():
    report = analyze("CCO", categories=["oxidation"])
    assert ids(report) == {"ox_pcc_primary", "ox_jones_primary"}
    assert ids(analyze("CCO", ids=["pbr3"])) == {"pbr3"}


@pytest.mark.parametrize("bad", ["", "C1CC", "not-a-molecule", "C(C)(C)(C)(C)C"])
def test_invalid_input(bad):
    with pytest.raises(InvalidMoleculeError):
        analyze(bad)


def test_parse_helpers():
    assert len(parse_mixture("CCO.CC(=O)O")) == 2
    assert [Chem.MolToSmiles(m) for m in parse_mixture(["CCO", "C"])] == ["CCO", "C"]
    assert Chem.MolToSmiles(parse_molecule(Chem.MolFromSmiles("OCC"))) == "CCO"
    with pytest.raises(InvalidMoleculeError):
        parse_molecule(42)
