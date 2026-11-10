"""Tour of the synth-path-helper Python API.

Run from the repository root after installing the package:

    python examples/demo.py
"""

from synth_path_helper import analyze, find_routes, predict_products

# 1. Analyse a molecule: functional groups and every applicable transformation
report = analyze("OCC(=O)c1ccccc1")  # 2-hydroxyacetophenone
for component in report.components:
    print(component.smiles, "->", [g.group.name for g in component.groups])

for pred in report.predictions[:5]:
    t = pred.transformation
    print(f"{t.id:<22} {', '.join(t.reagents):<45} -> {pred.main_product.smiles}")
    for caution in pred.cautions:
        print(f"{'':<24}[{caution.severity}] {caution.message}")

# 2. One transformation, one answer
print(predict_products("CC(C)=C", "hbr_markovnikov"))   # ['CC(C)(C)Br']
print(predict_products("CC(C)=C", "hydroboration"))     # ['CC(C)CO']

# 3. Two-component reaction (the mixture is given as dot-separated SMILES)
fischer = analyze("CC(=O)O.CCO").find("fischer_esterification")[0]
print(fischer.major[0].smiles)                           # 'CCOC(C)=O'

# 4. Multi-step route search
for route in find_routes("toluene", "O=C(O)c1ccc([N+](=O)[O-])cc1"):
    print(" -> ".join([route.start] + [s.product for s in route.steps]))
    print("    via", [s.transformation.id for s in route.steps])
