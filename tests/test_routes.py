from synth_path_helper import find_routes


def route_ids(route):
    return [s.transformation.id for s in route.steps]


def test_one_step_route():
    routes = find_routes("1-propanol", "CCC(=O)O")
    assert routes and len(routes[0]) == 1
    assert route_ids(routes[0]) == ["ox_jones_primary"]


def test_two_step_route_alkene_intermediate():
    routes = find_routes("cyclohexanol", "cyclohexene oxide")
    assert routes
    assert route_ids(routes[0]) == ["dehydration", "epoxidation"]
    assert routes[0].steps[0].product == "C1=CCCCC1"


def test_order_of_events_matters_on_arenes():
    # p-nitrobenzoic acid: nitrate first (CH3 directs para), then oxidise.
    # Oxidising first would give a meta director and the wrong isomer.
    routes = find_routes("toluene", "O=C(O)c1ccc([N+](=O)[O-])cc1")
    assert routes
    assert route_ids(routes[0]) == ["eas_nitration", "benzylic_oxidation"]


def test_trivial_and_impossible_routes():
    assert len(find_routes("CCO", "OCC")[0]) == 0
    assert find_routes("cyclohexanol", "cyclohexene oxide", max_steps=1) == []
    assert find_routes("CCO", "c1ccccc1", max_steps=2) == []


def test_two_component_step_with_partner():
    target = "CC(C)(O)c1ccccc1"
    assert find_routes("bromobenzene", target, max_steps=2) == []
    routes = find_routes("bromobenzene", target, max_steps=2, partners=["acetone"])
    assert routes and route_ids(routes[0]) == ["grignard_formation", "grignard_carbonyl"]
    assert routes[0].steps[1].co_reactants == ("CC(C)=O",)


def test_exclude_transformation():
    routes = find_routes("1-propanol", "CCC(=O)O", exclude=["ox_jones_primary"])
    assert routes and len(routes[0]) == 2
    assert "ox_jones_primary" not in route_ids(routes[0])


def test_route_to_dict():
    d = find_routes("CCO", "CC=O")[0].to_dict()
    assert d["steps"][0]["transformation"] == "ox_pcc_primary"
    assert d["steps"][0]["reagents"]
