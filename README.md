# synth-path-helper

A rule-based helper for organic chemistry. Give it a molecule (SMILES or a common name) and it finds the functional groups, predicts the products of the reactions those groups can undergo, and lists the reagents, catalysts and solvents to use.

```text
$ synth-path-helper ethanol          # output trimmed
...
  [ox_pcc_primary] PCC oxidation of a primary alcohol to an aldehyde
    Reagents       : PCC (pyridinium chlorochromate)
    Solvent(s)     : CH2Cl2 (anhydrous)
    Major product  : CC=O  [C2H4O, 44.05 g/mol]

  [ox_jones_primary] Jones oxidation of a primary alcohol to a carboxylic acid
    Reagents       : Jones reagent (CrO3 in aqueous H2SO4)
    Solvent(s)     : Acetone
    Major product  : CC(=O)O  [C2H4O2, 60.05 g/mol]
...
```

I started it during my first organic chemistry courses and labs at the Université de Montréal, to have the course's "functional group + reagent → product" rules in one place I could query. Technically it is SMARTS pattern matching on molecular graphs with [RDKit](https://www.rdkit.org/), plus a reaction library stored as plain JSON that you can extend.

## Project history

First version in March 2023. I wrote the code myself and asked Claude, an AI assistant, to review it. Until November 2026 it only existed on my laptop, without version control, which is why the repository starts with one commit. The README was rewritten with Claude's help when I published it.

## Features

- Functional-group perception: 29 groups defined as SMARTS patterns. It tells primary, secondary and tertiary alcohols apart, as well as alcohols from phenols, esters from ethers and amines from anilines.
- 69 named transformations across oxidation, reduction, additions, substitutions, eliminations, acid-derivative chemistry, Grignard/Wittig C–C bond formation and electrophilic aromatic substitution. Each one carries its reagents, catalysts, solvents, conditions, mechanism and alternatives.
- Major vs minor products, ranked by textbook selectivity rules: Markovnikov / anti-Markovnikov, Zaitsev / Hofmann and EAS directing effects (ortho/para vs meta, activation strength, steric hindrance).
- Chemoselectivity warnings. For example, "LiAlH4 will also reduce the ketone", "NaBH4 leaves the ester alone" or "Friedel–Crafts fails on a nitro-substituted ring".
- Excess reagent vs one equivalent: with several reactive sites (a diol, a diene) it shows both the fully converted product and the single-site products.
- Two-component reactions (Fischer esterification, Williamson ether synthesis, Grignard addition, reductive amination, and more). Pass a mixture such as `"CC(=O)O.CCO"`. If the partner is missing, the tool suggests one.
- Multi-step route search: a breadth-first search over the library finds short synthetic routes, e.g. toluene → *p*-nitrobenzoic acid (nitrate *first*, then oxidise).
- Common names for about 60 everyday compounds (`cyclohexanol`, `toluene`, `acetone`, …) in addition to SMILES.
- Text, JSON or SVG output; a Python API; user-extensible JSON libraries.

## Installation

Python 3.8–3.11 and RDKit 2022.03 or later. `requirements.txt` pins the versions I developed and tested with.

```bash
git clone https://github.com/jbenoitds/synth-path-helper.git
cd synth-path-helper
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt    # rdkit 2022.9.5, numpy 1.23.5
pip install -e .
```

This installs a `synth-path-helper` command. `python -m synth_path_helper` works too.

## Quick start (command line)

### Analyse a molecule

```bash
synth-path-helper analyze cyclohexanol          # 'analyze' is optional: synth-path-helper cyclohexanol
synth-path-helper analyze "CC(=O)c1ccccc1" --category Reduction
synth-path-helper analyze "OCC(=O)c1ccccc1" --brief --feasible
synth-path-helper analyze "CC(=O)O.CCO" --json  # a mixture: acid + alcohol
synth-path-helper analyze styrene --svg styrene.svg
```

Regiochemistry is decided by the rules, so isobutylene gives different products with HBr and with hydroboration:

```text
$ synth-path-helper analyze 2-methylpropene --only hbr_markovnikov --only hydroboration   # trimmed

  [hbr_markovnikov] Hydrobromination of an alkene (Markovnikov)
    Reagents       : HBr
    Solvent(s)     : CH2Cl2 or AcOH
    Major product  : CC(C)(C)Br  [C4H9Br, 137.02 g/mol]
    Minor          : CC(C)CBr  [C4H9Br, 137.02 g/mol]
    Selectivity    : Markovnikov: the new group ends up on the carbon that best
                     stabilises a positive charge.

  [hydroboration] Hydroboration-oxidation of an alkene (anti-Markovnikov, syn)
    Reagents       : 1) BH3·THF (or 9-BBN); 2) H2O2, NaOH(aq)
    Solvent(s)     : THF
    Major product  : CC(C)CO  [C4H10O, 74.12 g/mol]
    Minor          : CC(C)(C)O  [C4H10O, 74.12 g/mol]
```

| Option | Effect |
|---|---|
| `--category NAME` | only show one category (repeatable), e.g. `Oxidation` |
| `--only ID` | only run one transformation id (repeatable) |
| `--feasible` | hide reactions flagged as *not expected to work* |
| `--brief` | leave out conditions, mechanism and alternatives |
| `--json` | machine-readable output |
| `--svg FILE` | draw the input and the main product of every reaction |

### Find a synthetic route

```text
$ synth-path-helper route toluene "O=C(O)c1ccc([N+](=O)[O-])cc1"   # trimmed

Route 1 (2 steps)
  Step 1: Cc1ccccc1  ->  Cc1ccc([N+](=O)[O-])cc1
      Reaction       : Nitration of a benzene ring (EAS) [eas_nitration]
      Reagents       : conc. HNO3
      Catalyst(s)    : conc. H2SO4
  Step 2: Cc1ccc([N+](=O)[O-])cc1  ->  O=C(O)c1ccc([N+](=O)[O-])cc1
      Reaction       : Side-chain (benzylic) oxidation to a benzoic acid [benzylic_oxidation]
      Reagents       : KMnO4; NaOH (aq); then H3O+
```

Two-component steps need the partner molecule, which you supply with `--with`:

```bash
synth-path-helper route bromobenzene "CC(C)(O)c1ccccc1" --with acetone
#   Step 1: Brc1ccccc1 -> Br[Mg]c1ccccc1                     (Mg, Et2O)
#   Step 2: Br[Mg]c1ccccc1 + CC(C)=O -> CC(C)(O)c1ccccc1     (Grignard addition, H3O+)
```

Other route options: `--max-steps N` (default 3), `--max-routes N`, `--include-minor`, `--exclude ID`, `--json` and `--svg FILE`.

### Browse the library

```bash
synth-path-helper groups "OCC(=O)c1ccccc1"   # just the functional groups
synth-path-helper list                       # every transformation, by category
synth-path-helper show ox_pcc_primary        # the full record for one entry
```

## Quick start (Python)

```python
from synth_path_helper import analyze, predict_products, find_routes

report = analyze("CCOC(=O)CCC(C)=O")          # ethyl levulinate: ester + ketone
for comp in report.components:
    print([g.group.name for g in comp.groups])  # ['Ketone', 'Ester']

nabh4 = report.find("red_nabh4_ketone")[0]
print(nabh4.major[0].smiles)                  # CCOC(=O)CCC(C)O
print(nabh4.cautions[0].message)              # Ester is left untouched by NaBH4, ...

print(predict_products("CC(C)=C", "hbr_markovnikov"))   # ['CC(C)(C)Br']

for route in find_routes("cyclohexanol", "cyclohexene oxide"):
    print([s.transformation.id for s in route.steps])   # ['dehydration', 'epoxidation']

report.to_dict()                               # everything as JSON-ready dicts
```

`examples/demo.py` gives a longer tour.

## How it works

```
 SMILES / name ──► 1. parse & normalise ──► 2. perceive functional groups (SMARTS)
                                                      │
                                                      ▼
 report ◄── 6. format ◄── 5. cautions ◄── 4. rank ◄── 3. apply matching reaction templates
 (text / JSON / SVG)        (chemo-         (Markovnikov,   (RDKit RunReactants,
                            selectivity)     Zaitsev, EAS)   clean-up, de-duplication)
```

1. **Parsing** (`chem.py`). Input is looked up in the common-name dictionary, then parsed by RDKit and round-tripped through canonical SMILES. The round trip removes explicit-hydrogen flags (e.g. a typed `[CH2]`) that would otherwise carry wrong hydrogen counts into products. A dot-separated SMILES becomes several *components*.

2. **Functional-group perception** (`functional_groups.py`, `data/functional_groups.json`). Each group is a SMARTS pattern. For example, a secondary alcohol is

   ```
   [OX2H][CX4H1;!$(C([OX2H])[O,N,S,F,Cl,Br,I])]([#6])[#6]
   ```

   which reads: an OH on an sp3 CH that carries two carbons and no second heteroatom, so hemiacetals and hydrates do not count. Distinct matches are counted, so a diol shows "Primary alcohol ×2".

3. **Template matching** (`predictor.py`, `data/transformations.json`). Every transformation lists the groups that *trigger* it, plus one or more **reaction SMARTS** templates. PCC oxidation, for instance, is

   ```
   [CX4;H2,H3;!$(C([OX2H])[O,N,S,F,Cl,Br,I]):1][OX2H:2]>>[C:1]=[O:2]
   ```

   Templates are checked only against components that contain a trigger group. RDKit's `RunReactants` produces every possible product, one per matching site and orientation. Each product is cleaned up (hydrogens and charges recomputed on the atoms that changed, stale stereocentres cleared), sanitised and de-duplicated by canonical SMILES. Two-component templates are tried on every ordered pair of components. If a template has several sites and excess reagent makes sense, the template is applied repeatedly to give the fully converted product.

4. **Ranking** (`selectivity.py`). When a template can produce regioisomers, a named rule scores each product using the atom-map numbers that survive the reaction:

   | Rule | Score | Used by |
   |---|---|---|
   | `markovnikov` | how well the carbon that receives the new group stabilises a carbocation (alkyl substitution + benzylic bonus) | HBr, hydration, oxymercuration, bromohydrin, epoxide openings |
   | `anti_markovnikov` | the opposite | HBr/peroxides, hydroboration |
   | `zaitsev` | substitution of the new C=C (+ conjugation bonus) | acid dehydration, E2 with NaOEt |
   | `hofmann` | the opposite (less substituted wins) | E2 with t-BuOK |
   | `eas` | for each ring substituent: director type (o/p or meta) × strength × position, with an ortho penalty for bulky groups and a bonus at the naphthalene α-position | nitration, halogenation, sulfonation, Friedel–Crafts |

   Products scoring within 15% of the best one are reported together as major products. That is how a real ortho + para mixture (toluene) is shown, while *tert*-butylbenzene gives para only.

5. **Cautions**. Each entry can list other functional groups that cause trouble, each with a severity: `info` (chemoselective, good news), `warn` (that group reacts too) or `fails` (the reaction will not work, e.g. a Grignard reagent formed next to an OH). Predictions marked `fails` are shown with **NOT EXPECTED TO WORK** and are skipped by the route finder.

6. **Route search** (`routes.py`). This is a breadth-first search: every product (or product fragment) is a node, and every feasible prediction is an edge. The search stops at the first depth that reaches the target, compares molecules without stereochemistry, and lists the routes with the fewest warnings first.

## What the library covers

| Category | # | Transformations (ids) |
|---|---|---|
| Oxidation | 8 | `ox_pcc_primary`, `ox_jones_primary`, `ox_secondary`, `ox_aldehyde_acid`, `epoxidation`, `dihydroxylation`, `ozonolysis`, `benzylic_oxidation` |
| Reduction | 14 | `red_nabh4_aldehyde`, `red_nabh4_ketone`, `red_lialh4_ester`, `red_dibal_ester`, `red_lialh4_acid`, `red_lialh4_amide`, `red_lialh4_nitrile`, `red_wolff_kishner`, `red_nitro`, `hydrogenation_alkene`, `hydrogenation_alkyne`, `lindlar`, `dissolving_metal_alkyne`, `epoxide_lialh4` |
| Addition | 10 | `hbr_markovnikov`, `hbr_peroxide`, `hydration_acid`, `oxymercuration`, `hydroboration`, `bromination_alkene`, `bromohydrin`, `alkyne_hydration`, `alkyne_hydroboration`, `epoxide_acid_hydrolysis` |
| C–C bond formation | 5 | `acetylide_alkylation`, `sn2_cyanide`, `grignard_carbonyl`, `grignard_co2`, `wittig_methylenation` |
| Elimination | 3 | `dehydration`, `e2_ethoxide`, `e2_tbuok` |
| Substitution | 9 | `pbr3`, `socl2`, `hx_tertiary`, `williamson`, `sn2_hydroxide`, `sn2_azide`, `sn1_solvolysis`, `nbs_benzylic`, `epoxide_methoxide` |
| Protecting / activating groups | 3 | `tosylation`, `acetal_protection`, `amine_acetylation` |
| Organometallics | 1 | `grignard_formation` |
| Acid derivatives | 7 | `acid_to_acyl_chloride`, `fischer_esterification`, `saponification`, `acyl_chloride_alcohol`, `acyl_chloride_amine`, `amide_hydrolysis`, `nitrile_hydrolysis` |
| C–N bond formation | 1 | `reductive_amination` |
| Aromatic substitution | 8 | `eas_nitration`, `eas_bromination`, `eas_chlorination`, `eas_sulfonation`, `fc_acylation`, `sandmeyer_br`, `sandmeyer_cn`, `diazonium_hydrolysis` |

**Functional groups recognised:** primary / secondary / tertiary alcohols, phenols, aldehydes, ketones, carboxylic acids, esters, amides, acyl halides, anhydrides, nitriles, alkenes, alkynes (and terminal alkynes), benzene rings, primary / secondary / tertiary alkyl halides, aryl halides, Grignard reagents, primary / secondary / tertiary amines, anilines, nitro groups, ethers, epoxides and thiols.

Stereospecific outcomes are encoded where a template can express them: Lindlar gives the *cis* alkene and Na/NH3 the *trans* alkene. Syn/anti addition and inversion are described in each entry's conditions text.

## Adding your own reactions

Libraries are plain JSON. Load extra files with `--library` (repeatable) or `load_transformations([path])`. An entry with the same `id` as a built-in one replaces it. Every file is validated on load: SMARTS must parse, group ids must exist, and the number of templates must match the triggers.

```json
{
  "transformations": [
    {
      "id": "baeyer_villiger",
      "name": "Baeyer-Villiger oxidation of a cyclic ketone to a lactone",
      "category": "Oxidation",
      "triggers": ["ketone"],
      "smarts": "[#6:3][CX3;R:1](=[OX1:2])[#6;R:4]>>[#6:3][C:1](=[O:2])O[#6:4]",
      "reagents": ["mCPBA"],
      "catalysts": [],
      "solvents": ["CH2Cl2"],
      "conditions": "Room temperature.",
      "cautions": [
        {"if_present": ["alkene"], "message": "{group} would be epoxidised as well.", "severity": "warn"}
      ],
      "example": {"reactants": "O=C1CCCCC1", "major": "O=C1CCCCCO1"}
    }
  ]
}
```

```bash
synth-path-helper --library examples/custom_library.json analyze cyclohexanone
```

| Field | Meaning |
|---|---|
| `triggers` | group ids that make the entry relevant; for two-component reactions, a list of two lists (`[["carboxylic_acid"], ["alcohol_primary"]]`) |
| `smarts` | one reaction SMARTS, or a list of them. Several product fragments must be grouped as one product: `>>(A.B)` |
| `selectivity` | optional: `{"rule": "markovnikov", "site": 2, "other": 1}`, `{"rule": "zaitsev", "atoms": [1, 2]}`, `{"rule": "eas", "site": 1}` (numbers are atom maps) |
| `exhaustive` | `true` if excess reagent converts every matching site |
| `cautions` | `{"if_present": [...], "message": "... {group} ...", "severity": "info" \| "warn" \| "fails"}`; add `"always": true` to warn even when that group is one of the triggers |
| `partner_hint` | what to suggest when the second reactant is missing |
| `example` | `{"reactants": ..., "major": ...}`; the test suite checks that every entry predicts its own example |

## Project layout

```
synth-path-helper/
├── src/synth_path_helper/
│   ├── chem.py               # parsing, canonical SMILES, product clean-up
│   ├── functional_groups.py  # SMARTS-based group perception
│   ├── transformations.py    # loading/validating the reaction library
│   ├── selectivity.py        # Markovnikov, Zaitsev/Hofmann, EAS ranking
│   ├── predictor.py          # analyze(), predict(): the core engine
│   ├── routes.py             # breadth-first multi-step route search
│   ├── report.py             # text report and SVG drawing
│   ├── cli.py                # command-line interface
│   └── data/
│       ├── functional_groups.json
│       ├── transformations.json
│       └── common_compounds.json
├── tests/                    # pytest suite (249 tests)
├── examples/                 # demo.py, custom_library.json
├── .github/workflows/tests.yml
├── pyproject.toml
└── requirements*.txt
```

## Running the tests

```bash
pip install -r requirements-dev.txt
pip install -e .
python -m pytest -q
```

The suite has 249 tests. They cover group detection (including the near-misses: hemiacetals are not alcohols, esters are not ethers, anilines are not alkylamines), the documented example of every library entry, every selectivity rule, the cautions, two-component reactions, route search, the CLI and validation of malformed libraries. GitHub Actions runs it on Ubuntu and Windows with Python 3.8–3.11.

## Limitations

Keep in mind that this is a study aid, **not** a reaction predictor for the lab: the rules are the ones taught in undergraduate courses, and real chemistry has a lot more exceptions.

- No kinetics, thermodynamics, yields or real product ratios. "Major" means *favoured by the textbook rule*.
- Carbocation **rearrangements** (hydride or methyl shifts) are mentioned in the conditions but not generated.
- Stereochemistry is handled only where a template encodes it (Lindlar / Na-NH3). Other stereo outcomes (syn/anti, inversion, racemisation) are described in text, and stereocentres at the reaction centre are cleared.
- Cautions are checked per molecule, not per atom. A warning means the group is *present*, not that it sits next to the reaction site.
- SN2/E2 competition is simplified: substrate class picks the template, and a warning covers the borderline cases.
- The route finder explores only library reactions and stops at the shortest routes. It does not judge cost, protecting-group strategy beyond the few included, or practicality.

If you add a transformation, give it an `example`: the tests will check it automatically. Pull requests are welcome.

## License

MIT. See [LICENSE](LICENSE).
