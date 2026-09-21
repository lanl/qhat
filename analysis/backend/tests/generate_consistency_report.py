"""Generate comprehensive consistency test report."""

import json
import sys
from pathlib import Path


def load_results(filename):
    """Load test results from JSON."""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return []


def generate_report(results, output_file):
    """Generate markdown report from test results."""

    report = []
    report.append("# Backend Consistency Test Report")
    report.append("")
    report.append("This report summarizes consistency testing across QHAT backends.")
    report.append("")

    # Summary
    total_tests = len(results)
    passed_tests = sum(1 for r in results if r['data'].get('passed', False))

    report.append("## Summary")
    report.append("")
    report.append(f"- **Total tests**: {total_tests}")
    report.append(f"- **Passed**: {passed_tests}")
    report.append(f"- **Failed**: {total_tests - passed_tests}")
    report.append("")

    # Test details
    report.append("## Test Results")
    report.append("")

    for result in results:
        test_name = result['test_name']
        data = result['data']

        report.append(f"### {test_name}")
        report.append("")

        if 'backends_tested' in data:
            report.append(f"**Backends tested**: {', '.join(data['backends_tested'])}")
            report.append("")

        # Trotter tests
        if 'hamiltonian' in data:
            report.append(f"- **Hamiltonian**: {data['hamiltonian']}")
            report.append(f"- **Trotter order**: {data['trotter_order']}")
            report.append(f"- **Steps**: {data['num_steps']}")
            report.append(f"- **Qubits**: {data['num_qubits']}")
            report.append("")

            # Resources
            if 'resources' in data:
                report.append("**Resource Estimates:**")
                report.append("")
                report.append("| Backend | T Gates | Clifford | Rotations | Total |")
                report.append("|---------|---------|----------|-----------|-------|")
                for backend, res in data['resources'].items():
                    report.append(
                        f"| {backend} | {res['t_gates']} | {res['clifford_gates']} | "
                        f"{res['rotation_gates']} | {res['total_gates']} |"
                    )
                report.append("")

            # Matrix comparisons
            if 'matrix_comparisons' in data and data['matrix_comparisons']:
                report.append("**Matrix Consistency:**")
                report.append("")
                report.append("| Comparison | Rel. Diff. | Max Elem. Diff. | Status |")
                report.append("|------------|------------|-----------------|--------|")
                for comp in data['matrix_comparisons']:
                    rel_diff = comp['relative_difference']
                    max_diff = comp['max_element_difference']
                    status = "✅ Pass" if rel_diff < 1e-6 else "❌ Fail"
                    report.append(
                        f"| {comp['backends']} | {rel_diff:.2e} | {max_diff:.2e} | {status} |"
                    )
                report.append("")

        # Unitarity test
        elif 'unitarity_errors' in data:
            report.append("**Unitarity Errors (||U†U - I||):**")
            report.append("")
            report.append("| Backend | Error | Status |")
            report.append("|---------|-------|--------|")
            for backend, error in data['unitarity_errors'].items():
                status = "✅ Pass" if error < 1e-10 else "❌ Fail"
                report.append(f"| {backend} | {error:.2e} | {status} |")
            report.append("")

        # Controlled operations
        elif 'controlled_results' in data:
            report.append("**Controlled Operations:**")
            report.append("")
            report.append("| Backend | Qubits Increased | Identity Block Error | Status |")
            report.append("|---------|------------------|----------------------|--------|")
            for backend, res in data['controlled_results'].items():
                status = "✅ Pass" if res['identity_block_error'] < 1e-10 else "❌ Fail"
                report.append(
                    f"| {backend} | {res['qubits_increased']} | "
                    f"{res['identity_block_error']:.2e} | {status} |"
                )
            report.append("")

        # Adjoint operations
        elif 'adjoint_results' in data:
            report.append("**Adjoint Operations:**")
            report.append("")
            report.append("| Backend | ||U†U - I|| | ||adjoint() - U†|| | Status |")
            report.append("|---------|------------|---------------------|--------|")
            for backend, res in data['adjoint_results'].items():
                id_err = res['identity_error']
                adj_err = res['adjoint_correctness_error']
                status = "✅ Pass" if max(id_err, adj_err) < 1e-10 else "❌ Fail"
                report.append(
                    f"| {backend} | {id_err:.2e} | {adj_err:.2e} | {status} |"
                )
            report.append("")

        # Scaling test
        elif 'scaling_data' in data:
            report.append("**Resource Scaling with Step Count:**")
            report.append("")
            for backend, scaling in data['scaling_data'].items():
                if scaling:
                    report.append(f"*{backend}:*")
                    report.append("")
                    report.append("| Steps | T Gates | Total Gates |")
                    report.append("|-------|---------|-------------|")
                    for point in scaling:
                        report.append(
                            f"| {point['num_steps']} | {point['t_gates']} | "
                            f"{point['total_gates']} |"
                        )
                    report.append("")

        report.append(f"**Status**: {'✅ Passed' if data.get('passed', False) else '❌ Failed'}")
        report.append("")
        report.append("---")
        report.append("")

    # Write report
    with open(output_file, 'w') as f:
        f.write('\n'.join(report))

    print(f"Report generated: {output_file}")


if __name__ == "__main__":
    results_dir = Path("analysis/backend/tests/consistency_results")
    results_file = results_dir / "test_results.json"
    report_file = results_dir / "CONSISTENCY_REPORT.md"

    if not results_file.exists():
        print(f"Results file not found: {results_file}")
        print("Run the consistency tests first:")
        print("  pytest analysis/backend/tests/test_backend_consistency.py -v")
        sys.exit(1)

    results = load_results(results_file)

    if not results:
        print("No results found in results file")
        sys.exit(1)

    generate_report(results, report_file)
    print(f"\nGenerated report with {len(results)} test results")
