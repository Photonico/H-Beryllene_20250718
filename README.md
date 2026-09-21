# G-beryllene

Investigate the structural, electronic, optical, and topological properties of g-beryllene.

## Corrected optical results

The September 2026 optical audit fixes structure-derived thickness normalization,
absorption frequency units, and clipped figure ranges. The 34 requested figures,
per-file status, numerical evidence, and manuscript corrections are in
[exported_figures](exported_figures/README.md). The original VASP data are unchanged.

Install `requirements-optics.txt`, then run `python scripts/regenerate_optics.py`
and `python scripts/audit_optics.py` to regenerate and independently verify the
figures using the actual notebook code. Run `python -m pytest -q tests` for the
formula, geometry, plotting-entry and clipping regressions.

## License

The source code in the folder `vmatplot` is released under the MIT License.

You are welcome to use, modify, distribute, and adapt the code for academic, educational, or commercial purposes, provided that the original copyright and license notice are retained.

Unless otherwise stated, this license applies to the Python source code and related scripts in this repository. Research data, manuscript text, figures, and third-party files may be subject to separate terms.
