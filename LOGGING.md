# QHAT Logging Architecture

## Overview

QHAT uses a centralized logging architecture based on Python's standard `logging` module. This document explains how logging is configured and used throughout the codebase.

## Key Principles

1. **Configure Once, Use Everywhere**: Logging is configured once at application startup (in `driver.py` or `hamgen.py`)
2. **Module-Level Loggers**: Each module gets its own logger using `logging.getLogger(__name__)`
3. **Standard Python Patterns**: Uses Python's `logging` module idiomatically
4. **No Config Dependency**: Functions don't need `config_general` parameter just for logging

## Quick Start

### For Application Entry Points

```python
# In driver.py or hamgen.py
import logging
from qhat.common.logging_utils import configure_logging

# Load user configuration
state = load_configuration()

# Configure logging ONCE based on user settings
configure_logging(
    level=state.config_general.loglevel,
    logfile=state.config_general.logfile
)

# Now you can use standard logging
logger = logging.getLogger(__name__)
logger.info("Application started")
```

### For Module Code

```python
# In any module (e.g., analysis/unitary.py, common/trotter_flattened.py)
import logging

logger = logging.getLogger(__name__)

def my_function(param1, param2):
    logger.info(f"Starting computation with {param1}")
    logger.verbose(f"Detailed info: {param2}")
    logger.debug(f"Debug details: ...")
    
    # No need to pass config_general around!
```

## Log Levels

QHAT supports five log levels (in order of verbosity):

| Level | Numeric | Purpose | Example Use |
|-------|---------|---------|-------------|
| `DEBUG` | 10 | Detailed diagnostic information | Loop iterations, variable values |
| `VERBOSE` | 25 | Moderately detailed information | Intermediate computation results |
| `INFO` | 20 | General informational messages | "Starting computation", "Completed in 1.2s" |
| `WARNING` | 30 | Warning messages | Invalid config values, fallback behavior |
| `ERROR` | 40 | Error messages | Failed operations, exceptions |

### Custom VERBOSE Level

QHAT adds a custom `VERBOSE` level (25) between `INFO` (20) and `WARNING` (30) for moderately detailed logging:

```python
logger.verbose("This provides more detail than INFO but less than DEBUG")
```

## Configuration

### User Configuration

Users configure logging in their config files:

```python
# In analysis/config.py
general = GeneralConfigurationUser()
general.logfile = "my_analysis.log"
general.print_verbose()  # Sets level to "verbose"
# Or: general.print_debug() for "debug"
# Or: general.print_default() for "info"
```

### Programmatic Configuration

```python
from qhat.common.logging_utils import configure_logging

# Basic configuration
configure_logging(level="info", logfile="output.log")

# Custom format
configure_logging(
    level="verbose",
    logfile="output.log",
    format_string="{asctime} | {levelname} | {message}",
    include_module_name=False
)

# Console only (no file)
configure_logging(level="debug")
```

## Advanced Features

### Temporary Log Level Changes

Change log level for a specific section of code:

```python
from qhat.common.logging_utils import temporary_log_level
import logging

# Normal verbosity
logger.info("Normal logging")

# Temporarily enable debug logging
with temporary_log_level(logging.DEBUG):
    logger.debug("This is visible")
    some_detailed_computation()

# Back to normal
logger.debug("This is NOT visible")
```

### Execution Time Logging

Automatically log how long an operation takes:

```python
from qhat.common.logging_utils import log_execution_time

with log_execution_time("matrix computation", logger):
    result = expensive_computation()
# Automatically logs: "matrix computation completed in 1.23s"
```

### Dynamic Level Changes

Change log level during execution:

```python
from qhat.common.logging_utils import reconfigure_log_level

# Start at info level
configure_logging(level="info")

# Later, switch to verbose for debugging
reconfigure_log_level("verbose")
```

### Module-Specific Levels

Control verbosity per module:

```python
import logging

# Configure qhat logging at INFO level
configure_logging(level="info")

# But make trotter_coefficients_fast more verbose
logging.getLogger("qhat.analysis.trotter_coefficients_fast").setLevel(logging.DEBUG)
```

## Migration Guide

### Old Pattern (Before Refactoring)

```python
def encode_ramped_trotter(
        config_general: GeneralConfiguration,  # ❌ Just for logging!
        config_unitary: UnitaryConfiguration,
        hamiltonian,
        tevol_hbar):
    config_general.log("Starting encoding...")
    config_general.log_verbose("Detailed info...")
```

### New Pattern (After Refactoring)

```python
import logging

logger = logging.getLogger(__name__)

def encode_ramped_trotter(
        config_unitary: UnitaryConfiguration,  # ✅ Cleaner signature
        hamiltonian,
        tevol_hbar):
    logger.info("Starting encoding...")
    logger.verbose("Detailed info...")
```

### Backward Compatibility

During migration, some functions keep `config_general` as an optional parameter:

```python
def save_matrix(output_path, unitary_matrix, config_general=None, ...):
    """
    config_general: DEPRECATED - kept for backward compatibility
    """
    # Uses logger instead
    logger.info(f"Matrix saved to {output_path}")
```

## Best Practices

### DO ✅

- Create module-level logger: `logger = logging.getLogger(__name__)`
- Use appropriate log levels: `logger.info()` for key events, `logger.verbose()` for details
- Log operation completion with timing: `with log_execution_time(...)`
- Include context in messages: `logger.info(f"Computed {n_terms} terms in {time:.2f}s")`

### DON'T ❌

- Don't create loggers inside functions (use module-level)
- Don't pass `config_general` just for logging
- Don't use `print()` statements (use `logger.info()`)
- Don't concatenate strings with `+` (use f-strings: `f"value = {x}"`)
- Don't log sensitive information (credentials, keys, etc.)

## Log Output Format

### Default Format

```
2026-05-28 11:20:13,914    INFO | qhat.analysis.unitary          | Starting computation
2026-05-28 11:20:14,235 VERBOSE | qhat.analysis.unitary          | Intermediate result: 42
2026-05-28 11:20:15,891    INFO | qhat.analysis.unitary          | Computation completed in 1.98s
```

### Components

```
{timestamp:23} {level:>7s} | {module_name:35s} | {message}
```

- **Timestamp**: ISO format with milliseconds (23 chars)
- **Level**: Right-aligned, 7 chars (DEBUG, VERBOSE, INFO, WARNING, ERROR)
- **Module**: Fully-qualified Python module name (35 chars)
- **Message**: The log message

## Troubleshooting

### Logs Not Appearing

1. **Check log level**: Make sure your level is high enough
   ```python
   # DEBUG messages won't appear if level is INFO
   configure_logging(level="debug")  # Lower threshold
   ```

2. **Check logger name**: Make sure you're using a qhat logger
   ```python
   # Wrong: logging.getLogger("my_module")
   # Right: logging.getLogger("qhat.analysis.my_module")
   #   or: logging.getLogger(__name__)  # If in qhat package
   ```

3. **Check if logging is configured**: Call `configure_logging()` before logging
   ```python
   configure_logging(level="info")
   logger.info("Now this works")
   ```

### Too Much/Too Little Detail

Adjust the log level:

```python
# Too much detail? Raise the level
configure_logging(level="info")  # Hide verbose and debug

# Too little detail? Lower the level
configure_logging(level="verbose")  # Show more details
configure_logging(level="debug")    # Show everything
```

### File Not Created

Check that:
1. Directory exists or can be created
2. You have write permissions
3. Path is correct

```python
configure_logging(level="info", logfile="/tmp/test.log")  # Use absolute path
```

## Examples

### Complete Module Example

```python
"""
Example module demonstrating proper logging usage.
"""
import logging
from qhat.common.logging_utils import log_execution_time

logger = logging.getLogger(__name__)


def process_data(data):
    """Process some data with proper logging."""
    logger.info(f"Processing {len(data)} items")
    
    with log_execution_time("data processing", logger):
        results = []
        for i, item in enumerate(data):
            logger.verbose(f"Processing item {i+1}/{len(data)}")
            logger.debug(f"Item details: {item}")
            result = _process_item(item)
            results.append(result)
    
    logger.info(f"Successfully processed {len(results)} items")
    return results


def _process_item(item):
    """Helper function."""
    logger.debug(f"Processing: {item}")
    # ... processing logic ...
    return result
```

### Complete Application Example

```python
"""
Example application entry point.
"""
import logging
from qhat.common.logging_utils import configure_logging
from qhat.analysis.configuration import load_configuration

logger = logging.getLogger(__name__)


def main():
    # Load configuration
    state = load_configuration()
    
    # Configure logging
    configure_logging(
        level=state.config_general.loglevel,
        logfile=state.config_general.logfile
    )
    
    # Application startup
    logger.info("=" * 80)
    logger.info("QHAT Analysis Starting")
    logger.info("=" * 80)
    logger.info(f"Configuration: {state.config_script}")
    logger.info(f"Git hash: {state.config_general.git_hash}")
    
    try:
        # Run analysis
        result = run_analysis(state)
        logger.info("Analysis completed successfully")
        return result
    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
```

## API Reference

See `qhat/common/logging_utils.py` for complete API documentation.

### Main Functions

- `configure_logging(level, logfile, format_string, include_module_name)` - Configure logging
- `get_logger(name)` - Get a logger (prefer `logging.getLogger(__name__)`)
- `temporary_log_level(level, logger_name)` - Context manager for temporary level change
- `log_execution_time(operation, logger)` - Context manager for timing operations
- `reconfigure_log_level(level)` - Change log level dynamically

### Constants

- `VERBOSE` - Numeric value for VERBOSE log level (25)
