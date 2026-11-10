"""synth-path-helper: rule-based functional-group analysis and reaction prediction.

Quick start::

    >>> from synth_path_helper import analyze
    >>> report = analyze("CCO")                     # ethanol
    >>> report.find("ox_pcc_primary")[0].major[0].smiles
    'CC=O'
"""

__version__ = "0.1.0"

from .chem import InvalidMoleculeError, parse_molecule, parse_mixture  # noqa: E402
from .functional_groups import detect_groups, load_groups  # noqa: E402
from .predictor import Analysis, Prediction, Product, analyze, predict, predict_products  # noqa: E402
from .routes import Route, find_routes  # noqa: E402
from .transformations import Transformation, get_transformation, load_transformations  # noqa: E402

__all__ = [
    "__version__",
    "Analysis",
    "InvalidMoleculeError",
    "Prediction",
    "Product",
    "Route",
    "Transformation",
    "analyze",
    "detect_groups",
    "find_routes",
    "get_transformation",
    "load_groups",
    "load_transformations",
    "parse_mixture",
    "parse_molecule",
    "predict",
    "predict_products",
]
