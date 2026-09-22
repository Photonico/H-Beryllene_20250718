"""Portable fail-closed checks using tiny synthetic VASP files only."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
SCRIPT = Path(__file__).with_name("check_topology_outputs.py")
SPEC = importlib.util.spec_from_file_location("validator", SCRIPT)
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


class ValidatorTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="topology-validator-test-")
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def run_fixture(self, nelect=2, iteration=3, soc="T", footer=True,
                    truncated=False, atom_count=1, inversion_expected=True,
                    manifest_flat=False, eigen_nbands=6):
        run = self.root / "new-run"
        run.mkdir()
        (run / "POSCAR").write_text("Be slab\n1\n2 0 0\n0 2 0\n0 0 20\nBe\n1\nDirect\n0 0 .2\n")
        for name in ("INCAR", "KPOINTS", "POTCAR"):
            (run / name).write_text("Synthetic test input\n")
        (run / "OUTCAR").write_text(
            "vasp.6.5.0\nLSORBIT = " + soc + "\nLNONCOLLINEAR = T\nNELECT = " + str(nelect)
            + "\nNKPTS=4 NBANDS=6 NELM=10 NSW=0 NIONS=1\nIteration 1( " + str(iteration)
            + ")\naborting loop because EDIFF is reached\n"
            + ("General timing and accounting informations for this job:\n" if footer else ""))
        lines = ["1 1 1 1", "0", "0", "CAR", "synthetic", f"{nelect} 4 {eigen_nbands}"]
        for k in ((0, 0, 0), (.5, 0, 0), (0, .5, 0), (.5, .5, 0)):
            lines += ["", f"{k[0]} {k[1]} {k[2]} .25"]
            for band, energy in enumerate([-2, -2, 1, 1, 2, 2][:eigen_nbands], 1):
                lines.append(f"{band} {energy} {1 if band <= nelect else 0}")
        (run / "EIGENVAL").write_text("\n".join(lines[:-3] if truncated else lines) + "\n")
        entry = {"id": "alpha", "expected_nelect": nelect, "reference": "reference/path",
                 "atom_count": atom_count, "inversion_expected": inversion_expected}
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps(entry if manifest_flat else {"structures": [entry]}))
        args = [sys.executable, "-B", str(SCRIPT), str(run), "--manifest", str(manifest),
                "--structure-id", "alpha", "--no-spglib"]
        process = subprocess.run(args, capture_output=True, text=True)
        result = json.loads((run / "topology_validation.json").read_text())
        return process, result, args, run

    def test_completed_static_soc_and_no_overwrite(self):
        process, result, args, run = self.run_fixture()
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertTrue(result["eigenvalues"]["all_four_2d_trim_present"])
        self.assertEqual(result["eigenvalues"]["manifolds"][0]["direct_gap_ev"], 3)
        self.assertFalse(result["topology_certified"])
        path = run / "topology_validation.json"
        before = path.read_bytes()
        self.assertNotEqual(subprocess.run(args, capture_output=True).returncode, 0)
        self.assertEqual(before, path.read_bytes())

    def test_flat_manifest(self):
        process, result, _, _ = self.run_fixture(manifest_flat=True)
        self.assertEqual(process.returncode, 0, result)

    def test_odd_filling_keeps_physical_split_and_even_candidates_separate(self):
        process, result, _, _ = self.run_fixture(nelect=3)
        self.assertEqual(process.returncode, 0)
        manifolds = result["eigenvalues"]["manifolds"]
        self.assertEqual([m["N"] for m in manifolds], [2, 3, 4])
        self.assertFalse(manifolds[1]["time_reversal_compatible_dimension"])
        self.assertEqual(manifolds[1]["direct_gap_ev"], 0)

    def test_nelm_limit_skips_analysis(self):
        process, result, _, _ = self.run_fixture(iteration=10)
        self.assertEqual(process.returncode, 1)
        self.assertNotIn("eigenvalues", result)

    def test_missing_timing_footer(self):
        process, result, _, _ = self.run_fixture(footer=False)
        self.assertEqual(process.returncode, 1)
        self.assertFalse(result["validation_passed"])

    def test_non_soc(self):
        process, result, _, _ = self.run_fixture(soc="F")
        self.assertEqual(process.returncode, 1)
        self.assertNotIn("eigenvalues", result)

    def test_truncated_eigenval(self):
        process, result, _, _ = self.run_fixture(truncated=True)
        self.assertEqual(process.returncode, 1)
        self.assertTrue(any("Truncated" in error for error in result["errors"]))

    def test_manifest_atom_count_mismatch(self):
        process, result, _, _ = self.run_fixture(atom_count=2)
        self.assertEqual(process.returncode, 1)
        self.assertNotIn("eigenvalues", result)

    def test_manifest_inversion_mismatch(self):
        process, result, _, _ = self.run_fixture(inversion_expected=False)
        self.assertEqual(process.returncode, 1)
        self.assertNotIn("eigenvalues", result)

    def test_outcar_eigenval_dimensions_mismatch(self):
        process, result, _, _ = self.run_fixture(eigen_nbands=4)
        self.assertEqual(process.returncode, 1)
        self.assertTrue(any("NBANDS mismatch" in error for error in result["errors"]))

    def test_cartesian_and_negative_volume_poscar(self):
        path = self.root / "POSCAR"
        path.write_text("Be\n-8\n1 0 0\n0 1 0\n0 0 1\nBe\n1\nCartesian\n.25 .5 .75\n")
        structure = VALIDATOR.read_poscar(path)
        self.assertEqual(structure["lattice"], [[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]])
        self.assertEqual(structure["fractional_positions"], [[.25, .5, .75]])

    def test_inverse_and_skew_minimum_image(self):
        lattice = [[2., 0., 0.], [1.9, .2, 0.], [0., 0., 20.]]
        inverse = VALIDATOR.inverse(lattice)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(sum(lattice[i][k] * inverse[k][j] for k in range(3)), float(i == j))
        # Componentwise fractional rounding gives a much longer image here.
        distance = VALIDATOR.periodic_distance([.49, .49, 0.], lattice, inverse)
        brute = min(sum(v * v for v in VALIDATOR.row_product([.49 + i, .49 + j, 0.], lattice)) ** .5
                    for i in range(-5, 6) for j in range(-5, 6))
        self.assertAlmostEqual(distance, brute)


if __name__ == "__main__":
    unittest.main()
