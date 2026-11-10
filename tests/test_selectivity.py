import pytest

from synth_path_helper import analyze, predict_products
from synth_path_helper.chem import canonicalize
from synth_path_helper.selectivity import is_major


def majors(smiles, tid):
    return set(predict_products(smiles, tid))


def minors(smiles, tid):
    pred = analyze(smiles, ids=[tid]).find(tid)[0]
    return {p.smiles for p in pred.minor}


def c(s):
    return canonicalize(s)


# --- addition to alkenes ---------------------------------------------------

def test_markovnikov_tertiary_carbon():
    assert majors("CC(C)=C", "hbr_markovnikov") == {c("CC(C)(C)Br")}
    assert minors("CC(C)=C", "hbr_markovnikov") == {c("CC(C)CBr")}


def test_markovnikov_benzylic_carbon():
    assert majors("C=Cc1ccccc1", "hbr_markovnikov") == {c("CC(Br)c1ccccc1")}


def test_markovnikov_no_preference_gives_mixture():
    # Both carbons of 2-pentene are secondary, so both regioisomers are major.
    assert majors("CC=CCC", "hbr_markovnikov") == {c("CCC(Br)CC"), c("CCCC(C)Br")}


def test_anti_markovnikov_hydroboration():
    assert majors("CC1=CCCCC1", "hydroboration") == {c("CC1CCCCC1O")}
    assert majors("C=Cc1ccccc1", "hydroboration") == {c("OCCc1ccccc1")}


def test_bromohydrin_oh_on_more_substituted_carbon():
    assert majors("CC(C)=C", "bromohydrin") == {c("CC(C)(O)CBr")}


# --- eliminations ----------------------------------------------------------

@pytest.mark.parametrize(
    "smiles, tid, expected",
    [
        ("CCC(C)(C)Br", "e2_ethoxide", "CC=C(C)C"),     # Zaitsev: trisubstituted
        ("CCC(C)(C)Br", "e2_tbuok", "C=C(C)CC"),        # Hofmann: disubstituted
        ("CC(Br)CC", "e2_ethoxide", "CC=CC"),
        ("CC(Br)CC", "e2_tbuok", "C=CCC"),
        ("OC1CCCCC1C", "dehydration", "CC1=CCCCC1"),
        ("CC(O)Cc1ccccc1", "dehydration", "CC=Cc1ccccc1"),  # conjugation with the ring wins
    ],
)
def test_elimination_regiochemistry(smiles, tid, expected):
    assert majors(smiles, tid) == {c(expected)}


# --- electrophilic aromatic substitution -----------------------------------

def test_meta_director():
    assert majors("O=[N+]([O-])c1ccccc1", "eas_nitration") == {c("O=[N+]([O-])c1cccc([N+](=O)[O-])c1")}
    assert majors("CC(=O)c1ccccc1", "eas_bromination") == {c("CC(=O)c1cccc(Br)c1")}


def test_ortho_para_director_gives_o_p_mixture():
    found = majors("Cc1ccccc1", "eas_nitration")
    assert c("Cc1ccc([N+](=O)[O-])cc1") in found       # para
    assert c("Cc1ccccc1[N+](=O)[O-]") in found         # ortho
    assert c("Cc1cccc([N+](=O)[O-])c1") in minors("Cc1ccccc1", "eas_nitration")  # meta is minor


def test_bulky_group_favours_para():
    assert majors("CC(C)(C)c1ccccc1", "eas_nitration") == {c("CC(C)(C)c1ccc([N+](=O)[O-])cc1")}


def test_strongest_activator_wins():
    # p-nitrotoluene: CH3 (o/p) and NO2 (meta) both point at C2
    assert majors("Cc1ccc([N+](=O)[O-])cc1", "eas_nitration") == {c("Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]")}
    # p-methylanisole: OMe beats CH3, so substitution goes ortho to OMe
    assert majors("COc1ccc(C)cc1", "eas_bromination") == {c("COc1ccc(C)cc1Br")}


def test_naphthalene_alpha_position():
    assert majors("c1ccc2ccccc2c1", "eas_bromination") == {c("Brc1cccc2ccccc12")}


def test_halogens_are_ortho_para_directors():
    found = majors("Clc1ccccc1", "eas_nitration")
    assert c("O=[N+]([O-])c1ccc(Cl)cc1") in found
    assert c("O=[N+]([O-])c1cccc(Cl)c1") not in found


def test_unsubstituted_benzene_single_product():
    assert majors("c1ccccc1", "eas_nitration") == {c("O=[N+]([O-])c1ccccc1")}


# --- epoxide opening ---------------------------------------------------------

def test_basic_epoxide_opening_at_less_hindered_carbon():
    assert majors("CC1(C)CO1", "epoxide_methoxide") == {c("COCC(C)(C)O")}
    assert majors("CC1(C)CO1", "epoxide_lialh4") == {c("CC(C)(C)O")}


def test_is_major_tolerance():
    assert is_major(2.0, 2.0)
    assert is_major(1.75, 2.0)     # ortho vs para for toluene
    assert not is_major(0.6, 2.0)  # ortho vs para for tert-butylbenzene
    assert is_major(0.0, 0.0)
    assert not is_major(-2.0, 0.0)
