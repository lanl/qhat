#!/usr/bin/env python3.11
"""
Script to check the status of logging migration.

This script scans the codebase and reports:
1. Functions that still use config_general for logging
2. Files that need to import logging module
3. Usage statistics

Run from the qhat root directory:
    python3.11 scripts/check_logging_migration.py
"""

import os
import re
from pathlib import Path
from collections import defaultdict


def find_python_files(root_dir):
    """Find all Python files in the project."""
    python_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Skip certain directories
        if any(skip in dirpath for skip in ['__pycache__', '.git', 'venv', 'env']):
            continue
        for filename in filenames:
            if filename.endswith('.py'):
                python_files.append(Path(dirpath) / filename)
    return python_files


def analyze_file(filepath):
    """Analyze a Python file for logging patterns."""
    with open(filepath, 'r') as f:
        content = f.read()

    results = {
        'has_logging_import': False,
        'has_logger_creation': False,
        'config_general_param_count': 0,
        'config_general_log_calls': 0,
        'standard_log_calls': 0,
        'functions_with_config_general': []
    }

    # Check for logging import
    if re.search(r'^import logging', content, re.MULTILINE):
        results['has_logging_import'] = True

    # Check for logger creation
    if re.search(r'logger\s*=\s*logging\.getLogger', content):
        results['has_logger_creation'] = True

    # Find functions with config_general parameter
    func_pattern = r'def\s+(\w+)\s*\([^)]*config_general[^)]*\):'
    for match in re.finditer(func_pattern, content):
        results['config_general_param_count'] += 1
        results['functions_with_config_general'].append(match.group(1))

    # Count config_general.log* calls
    log_call_pattern = r'config_general\.log[_\w]*\s*\('
    results['config_general_log_calls'] = len(re.findall(log_call_pattern, content))

    # Count standard logger.* calls
    standard_pattern = r'logger\.(debug|verbose|info|warning|error|critical)\s*\('
    results['standard_log_calls'] = len(re.findall(standard_pattern, content))

    return results


def print_summary(stats):
    """Print summary statistics."""
    print("\n" + "=" * 80)
    print("LOGGING MIGRATION STATUS")
    print("=" * 80)

    print(f"\nFiles analyzed: {stats['total_files']}")
    print(f"Files with config_general logging: {stats['files_with_old_logging']}")
    print(f"Files with standard logging: {stats['files_with_new_logging']}")
    print(f"Files fully migrated: {stats['files_migrated']}")

    migration_pct = (stats['files_migrated'] / stats['total_files'] * 100) if stats['total_files'] > 0 else 0
    print(f"\nMigration progress: {migration_pct:.1f}%")

    print(f"\nTotal config_general parameters: {stats['total_config_general_params']}")
    print(f"Total config_general.log* calls: {stats['total_old_log_calls']}")
    print(f"Total logger.* calls: {stats['total_new_log_calls']}")


def print_details(file_results):
    """Print detailed results for files needing migration."""
    needs_migration = [
        (filepath, results) for filepath, results in file_results.items()
        if results['config_general_log_calls'] > 0 or results['config_general_param_count'] > 0
    ]

    if not needs_migration:
        print("\n✅ All files have been migrated!")
        return

    print("\n" + "-" * 80)
    print("FILES NEEDING MIGRATION")
    print("-" * 80)

    for filepath, results in sorted(needs_migration):
        rel_path = filepath.relative_to(Path.cwd())
        print(f"\n📄 {rel_path}")

        if results['config_general_param_count'] > 0:
            print(f"   ⚠️  {results['config_general_param_count']} functions with config_general parameter:")
            for func_name in results['functions_with_config_general']:
                print(f"      - {func_name}()")

        if results['config_general_log_calls'] > 0:
            print(f"   ⚠️  {results['config_general_log_calls']} config_general.log* calls")

        if not results['has_logging_import']:
            print(f"   ℹ️  Needs: import logging")

        if not results['has_logger_creation']:
            print(f"   ℹ️  Needs: logger = logging.getLogger(__name__)")


def main():
    """Main function."""
    print("Scanning codebase for logging patterns...")

    # Find all Python files
    root_dir = Path.cwd()
    if not (root_dir / 'analysis').exists():
        print("Error: Run this script from the qhat root directory")
        return 1

    python_files = find_python_files(root_dir)
    print(f"Found {len(python_files)} Python files")

    # Analyze each file
    file_results = {}
    for filepath in python_files:
        # Skip test and example files for now
        if 'test_logging' in str(filepath) or 'check_logging_migration' in str(filepath):
            continue
        file_results[filepath] = analyze_file(filepath)

    # Compute statistics
    stats = {
        'total_files': len(file_results),
        'files_with_old_logging': 0,
        'files_with_new_logging': 0,
        'files_migrated': 0,
        'total_config_general_params': 0,
        'total_old_log_calls': 0,
        'total_new_log_calls': 0,
    }

    for results in file_results.values():
        if results['config_general_log_calls'] > 0 or results['config_general_param_count'] > 0:
            stats['files_with_old_logging'] += 1
        if results['standard_log_calls'] > 0:
            stats['files_with_new_logging'] += 1
        if results['standard_log_calls'] > 0 and results['config_general_log_calls'] == 0:
            stats['files_migrated'] += 1

        stats['total_config_general_params'] += results['config_general_param_count']
        stats['total_old_log_calls'] += results['config_general_log_calls']
        stats['total_new_log_calls'] += results['standard_log_calls']

    # Print results
    print_summary(stats)
    print_details(file_results)

    print("\n" + "=" * 80)
    print("For migration instructions, see LOGGING.md")
    print("=" * 80 + "\n")

    return 0


if __name__ == "__main__":
    exit(main())
