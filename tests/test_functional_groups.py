import pytest

from synth_path_helper import detect_groups, load_groups, parse_molecule
from synth_path_helper.functional_groups import get_group, group_ids


def ids(smiles):
    return group_ids(parse_molecule(smiles))


@pytest.mark.parametrize(
    "smiles, present, absent",
    [
        ("CCO", {"alcohol_primary"}, {"alcohol_secondary", "alcohol_tertiary", "phenol", "ether"}),
        ("CO", {"alcohol_primary"}, set()),
        ("CC(C)O", {"alcohol_secondary"}, {"alcohol_primary", "alcohol_tertiary"}),
        ("CC(C)(C)O", {"alcohol_tertiary"}, {"alcohol_primary", "alcohol_secondary"}),
        ("Oc1ccccc1", {"phenol", "arene"}, {"alcohol_primary", "alcohol_secondary"}),
        ("CC(O)O", set(), {"alcohol_secondary"}),  # gem-diol (aldehyde hydrate) is not an alcohol here
        ("CC=O", {"aldehyde"}, {"ketone"}),
        ("C=O", {"aldehyde"}, {"ketone", "alkene"}),
        ("CC(C)=O", {"ketone"}, {"aldehyde"}),
        ("CC(=O)O", {"carboxylic_acid"}, {"aldehyde", "ketone", "ester", "alcohol_primary"}),
        ("O=CO", {"carboxylic_acid"}, {"aldehyde"}),
        ("CCOC(C)=O", {"ester"}, {"ether", "ketone", "carboxylic_acid"}),
        ("CC(=O)OC(C)=O", {"anhydride"}, {"ester", "ether"}),
        ("CC(N)=O", {"amide"}, {"amine_primary", "ketone"}),
        ("CC(=O)Cl", {"acyl_halide"}, {"ketone", "alkyl_halide_primary"}),
        ("CC#N", {"nitrile"}, {"alkyne"}),
        ("C=CC", {"alkene"}, {"arene"}),
        ("c1ccccc1", {"arene"}, {"alkene"}),
        ("C#CC", {"alkyne", "alkyne_terminal"}, set()),
        ("CC#CC", {"alkyne"}, {"alkyne_terminal"}),
        ("CCBr", {"alkyl_halide_primary"}, {"alkyl_halide_secondary", "alkyl_halide_tertiary"}),
        ("CC(C)Br", {"alkyl_halide_secondary"}, {"alkyl_halide_primary"}),
        ("CC(C)(C)Br", {"alkyl_halide_tertiary"}, {"alkyl_halide_secondary"}),
        ("Brc1ccccc1", {"aryl_halide"}, {"alkyl_halide_primary"}),
        ("CCN", {"amine_primary"}, {"amine_secondary", "aniline"}),
        ("CNC", {"amine_secondary"}, {"amine_primary"}),
        ("CN(C)C", {"amine_tertiary"}, {"amine_secondary"}),
        ("Nc1ccccc1", {"aniline"}, {"amine_primary"}),
        ("CC(=O)Nc1ccccc1", {"amide"}, {"aniline", "ketone"}),
        ("O=[N+]([O-])c1ccccc1", {"nitro", "arene"}, set()),
        ("CCOCC", {"ether"}, {"epoxide"}),
        ("C1CO1", {"epoxide"}, {"ether"}),
        ("CCS", {"thiol"}, set()),
        ("CC[Mg]Br", {"grignard"}, {"alkyl_halide_primary"}),
    ],
)
def test_group_detection(smiles, present, absent):
    found = ids(smiles)
    assert present <= found, f"{smiles}: missing {present - found}"
    assert not (absent & found), f"{smiles}: unexpected {absent & found}"


def test_counts_distinct_occurrences():
    groups = {g.group.id: g for g in detect_groups(parse_molecule("OCCO"))}
    assert groups["alcohol_primary"].count == 2
    naphthalene = {g.group.id: g for g in detect_groups(parse_molecule("c1ccc2ccccc2c1"))}
    assert naphthalene["arene"].count == 2


def test_multifunctional_molecule():
    # 4-hydroxy-4-phenylbutan-2-one: ketone + secondary alcohol + benzene ring
    assert {"ketone", "alcohol_secondary", "arene"} <= ids("CC(=O)CC(O)c1ccccc1")


def test_library_is_well_formed():
    groups = load_groups()
    assert len({g.id for g in groups}) == len(groups)
    for g in groups:
        assert g.pattern is not None
        assert g.name and g.family and g.notes
    assert get_group("ketone").name == "Ketone"
    with pytest.raises(KeyError):
        get_group("not_a_group")


def test_to_dict():
    d = detect_groups(parse_molecule("CCO"))[0].to_dict()
    assert d["id"] == "alcohol_primary" and d["count"] == 1 and d["atoms"]
