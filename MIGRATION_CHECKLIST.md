# Logging Migration Checklist

Quick reference for migrating files to the new logging architecture.

## Pre-Migration

- [ ] Read LOGGING.md (at least the Quick Start section)
- [ ] Run migration checker: `python3.11 scripts/check_logging_migration.py`
- [ ] Pick a file to migrate (start with small files)
- [ ] Make sure you're on the `refactor/logging-architecture` branch

## For Each File

### 1. Add Imports
```python
import logging

logger = logging.getLogger(__name__)
```

- [ ] Add `import logging` at top of file
- [ ] Add `logger = logging.getLogger(__name__)` after imports

### 2. Replace Log Calls

Replace each occurrence:

```python
# Before → After
config_general.log()          → logger.info()
config_general.log_verbose()  → logger.verbose()
config_general.log_debug()    → logger.debug()
```

- [ ] Replace all `config_general.log()` → `logger.info()`
- [ ] Replace all `config_general.log_verbose()` → `logger.verbose()`
- [ ] Replace all `config_general.log_debug()` → `logger.debug()`
- [ ] Check for warnings: `config_general.log("WARNING: ...")` → `logger.warning(...)`

### 3. Update Function Signatures

**Only if** `config_general` is ONLY used for logging:

```python
# Before
def my_function(config_general: GeneralConfiguration, other_params):
    config_general.log("message")
    return result

# After
def my_function(other_params):
    logger.info("message")
    return result
```

⚠️ **Don't remove if used for other things**:
- `config_general.git_hash`
- `config_general.t_hbar`
- Any other non-logging attributes

If config_general is used for both logging and other things:
- [ ] Keep the parameter
- [ ] Just replace the log calls
- [ ] Add TODO comment: `# TODO: Remove config_general once all uses eliminated`

### 4. Update Callers

If you removed `config_general` parameter:

```python
# Before
result = my_function(config_general, other_args)

# After
result = my_function(other_args)
```

- [ ] Find all callers of modified functions
- [ ] Remove `config_general` argument from calls
- [ ] Make sure other arguments are in correct position

### 5. Test

- [ ] Run the module: `python3.11 analysis/driver.py` (or relevant script)
- [ ] Check log output appears correctly
- [ ] Verify log levels work (try different levels)
- [ ] Run relevant tests: `pytest analysis/tests/test_*.py`

### 6. Verify

- [ ] Run migration checker again: `python3.11 scripts/check_logging_migration.py`
- [ ] Confirm this file no longer listed (or fewer issues)
- [ ] Check diff: `git diff analysis/your_file.py`
- [ ] Make sure changes look correct

## Commit

```bash
git add analysis/your_file.py
git commit -m "Migrate your_file.py to new logging architecture

- Replace config_general.log* with logger.*
- Remove config_general parameter from function_name()
- Update N callers

Remaining in this module: 0 config_general.log calls"
```

- [ ] Stage file: `git add <file>`
- [ ] Write clear commit message
- [ ] Mention number of calls migrated
- [ ] Note any remaining issues

## Common Patterns

### Pattern 1: Simple Replacement

```python
# Before
def compute(config_general, data):
    config_general.log(f"Processing {len(data)} items")
    result = process(data)
    config_general.log_verbose(f"Result: {result}")
    return result

# After
import logging
logger = logging.getLogger(__name__)

def compute(data):
    logger.info(f"Processing {len(data)} items")
    result = process(data)
    logger.verbose(f"Result: {result}")
    return result
```

### Pattern 2: config_general Used for Other Things

```python
# Before
def encode(config_general, config_unitary, hamiltonian):
    config_general.log("Starting encoding")
    timestep = config_general.t_hbar  # ⚠️ Non-logging use!
    # ...

# After - Keep parameter!
import logging
logger = logging.getLogger(__name__)

def encode(config_general, config_unitary, hamiltonian):
    logger.info("Starting encoding")
    timestep = config_general.t_hbar  # Still needed
    # ...
```

### Pattern 3: Optional config_general for Backward Compatibility

```python
# Before
def save_matrix(path, matrix, config_general):
    # ... save logic ...
    config_general.log(f"Saved to {path}")

# After
import logging
logger = logging.getLogger(__name__)

def save_matrix(path, matrix, config_general=None):
    """
    config_general: DEPRECATED - kept for backward compatibility
    """
    # ... save logic ...
    logger.info(f"Saved to {path}")
```

### Pattern 4: Warnings

```python
# Before
config_general.log(f"WARNING: Invalid value '{value}', using default")

# After  
logger.warning(f"Invalid value '{value}', using default")
```

## Priority Order

Migrate in this order for easiest progression:

1. **Small utility files** (easiest)
   - file_io.py ✅ Done
   - ordering.py
   
2. **Core analysis files**
   - algorithm.py (9 calls)
   - analysis.py (12 calls)
   - hamiltonian.py (19 calls)
   
3. **Large compute files** (hardest)
   - trotter_coefficients_fast.py (44 calls)

4. **Common modules**
   - common/*.py files

5. **Hamiltonian generator**
   - hamiltonian_generator/*.py files

6. **Infrastructure cleanup**
   - config_types.py (remove logging methods)
   - Remove old _configure_log() function

## Troubleshooting

### "NameError: name 'logger' is not defined"
- [ ] Did you add `logger = logging.getLogger(__name__)`?
- [ ] Is it at module level (not inside a function)?

### "AttributeError: 'Logger' object has no attribute 'verbose'"
- [ ] Did you call `configure_logging()` at startup?
- [ ] Is it called before any logging happens?
- [ ] Check driver.py has the configure_logging() call

### "No logs appearing"
- [ ] Is logging configured? Check driver.py
- [ ] Is log level too high? Try `level="debug"`
- [ ] Are you using correct logger name? Should be under "qhat.*"

### "Tests failing"
- [ ] Did you update test fixtures?
- [ ] Do tests need to configure logging?
- [ ] Did you update function calls to match new signature?

## Quick Reference Commands

```bash
# Check migration status
python3.11 scripts/check_logging_migration.py

# Find all config_general.log calls in a file
grep -n "config_general\.log" analysis/your_file.py

# Find all functions with config_general parameter
grep -n "def.*config_general" analysis/your_file.py

# Test logging works
python3.11 test_logging.py

# Run migration checker for specific file
python3.11 scripts/check_logging_migration.py | grep "your_file.py" -A 10
```

## Resources

- **LOGGING.md** - Complete documentation
- **REFACTOR_SUMMARY.md** - Overall status and plan
- **logging_utils.py** - API reference
- **test_logging.py** - Working examples
- **scripts/check_logging_migration.py** - Track progress

## Getting Help

If stuck, check:
1. LOGGING.md "Troubleshooting" section
2. LOGGING.md "Migration Guide" section
3. Look at already-migrated files:
   - analysis/unitary.py (complete example)
   - analysis/file_io.py (simple example)
4. Run the migration checker for guidance
5. Look at test_logging.py for examples

---

**Remember**: This is a gradual migration. It's okay to keep `config_general` parameter if it's used for non-logging purposes. Just replace the log calls first!
