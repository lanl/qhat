import tempfile
import unittest
import zipfile
from pathlib import Path

from first_reproducibility_experiment import compare_metric, copy_verified, digest, json_digest, save, validate_tensor_archive


class AuditTests(unittest.TestCase):
    def test_historical_f2_difference_is_flagged(self):
        self.assertEqual(compare_metric("one_minus_overlap", 5.8586e-11, 1.0582e-7,
                                        repeat=False), "different")

    def test_counts_are_exact(self):
        self.assertEqual(compare_metric("number_of_pauli_terms", 5064, 5098,
                                        repeat=True), "different")

    def test_small_noise_allowed(self):
        self.assertEqual(compare_metric("one_minus_overlap", 2.8558e-10, 2.8559e-10,
                                        repeat=True), "match")

    def test_floor_is_explicit(self):
        self.assertEqual(compare_metric("one_minus_overlap", 0., 1e-14,
                                        repeat=False), "below_floor_match")

    def test_nonfinite_rejected(self):
        self.assertEqual(compare_metric("bch2_hf_state_norm", float("nan"), .1,
                                        repeat=True), "nonfinite")

    def test_missing_constant_rejected_in_preflight(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.tensors.npz"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("one_body.npy", b"")
                archive.writestr("two_body.npy", b"")
            with self.assertRaisesRegex(ValueError, "constant.npy"):
                validate_tensor_archive(path)

    def test_frozen_copy_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, target = root / "input", root / "frozen" / "input"
            source.write_bytes(b"unchanged original")
            self.assertEqual(copy_verified(source, target), digest(source))
            self.assertEqual(source.read_bytes(), target.read_bytes())
            with self.assertRaises(FileExistsError):
                copy_verified(source, target)

    def test_result_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            save(path, {"a": 1})
            with self.assertRaises(FileExistsError):
                save(path, {"a": 2})

    def test_fingerprint_is_order_and_coefficient_sensitive(self):
        self.assertNotEqual(json_digest([["X", .1], ["Z", .2]]),
                            json_digest([["Z", .2], ["X", .1]]))
        self.assertNotEqual(json_digest([["X", .1]]), json_digest([["X", .2]]))


if __name__ == "__main__":
    unittest.main()
