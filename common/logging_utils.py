"""
Centralized logging utilities for QHAT.

This module provides a standardized logging setup for the QHAT toolkit.
Configure logging once at application startup, then use standard module-level
loggers throughout the codebase.

Usage:
    # In main entry point (driver.py, hamgen.py):
    from qhat.logging_utils import configure_logging
    configure_logging(level="verbose", logfile="analysis.log")

    # In any module:
    import logging
    logger = logging.getLogger(__name__)
    logger.info("Starting computation...")
    logger.verbose("Detailed information...")
"""

import logging
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Optional


# Custom log levels
# Lower values = more verbose in Python logging
VERBOSE = 15  # Between DEBUG (10) and INFO (20)
RESULTS = 45  # Between ERROR (40) and CRITICAL (50) - for critical output that should almost always be shown


def add_custom_levels():
    """
    Add custom log levels to logging module.

    VERBOSE: Between DEBUG and INFO for moderately detailed logging.
    RESULTS: Between ERROR and CRITICAL for critical output (results, completion)
             that should be shown unless user explicitly sets level to CRITICAL.
    """
    # Add VERBOSE level
    if not hasattr(logging, 'VERBOSE'):
        logging.addLevelName(VERBOSE, "VERBOSE")

        def verbose(self, message, *args, **kwargs):
            if self.isEnabledFor(VERBOSE):
                self._log(VERBOSE, message, args, **kwargs)

        logging.Logger.verbose = verbose
        logging.verbose = lambda msg, *args, **kwargs: logging.log(VERBOSE, msg, *args, **kwargs)
        logging.VERBOSE = VERBOSE

    # Add RESULTS level
    if not hasattr(logging, 'RESULTS'):
        logging.addLevelName(RESULTS, "RESULTS")

        def results(self, message, *args, **kwargs):
            if self.isEnabledFor(RESULTS):
                self._log(RESULTS, message, args, **kwargs)

        logging.Logger.results = results
        logging.results = lambda msg, *args, **kwargs: logging.log(RESULTS, msg, *args, **kwargs)
        logging.RESULTS = RESULTS


def configure_logging(
    level: str = "info",
    logfile: Optional[str] = None,
    format_string: Optional[str] = None,
    include_module_name: bool = True
) -> logging.Logger:
    """
    Configure logging for QHAT application.

    Call this ONCE at application startup (e.g., in driver.py or hamgen.py).
    After calling this, all modules can use standard logging:
        logger = logging.getLogger(__name__)
        logger.info("message")

    Args:
        level: Log level string - one of: "debug", "verbose", "info", "warning", "error"
        logfile: Path to log file (if None, only logs to console)
        format_string: Custom format string (uses default if None)
        include_module_name: If True, includes module name in log messages

    Returns:
        The configured qhat root logger

    Raises:
        ValueError: If log level is invalid

    Example:
        >>> from qhat.logging_utils import configure_logging
        >>> configure_logging(level="verbose", logfile="analysis.log")
        >>> import logging
        >>> logger = logging.getLogger(__name__)
        >>> logger.info("Application started")
    """
    # Add custom levels (VERBOSE and RESULTS)
    add_custom_levels()

    # Parse level string
    level_map = {
        "debug": logging.DEBUG,
        "verbose": VERBOSE,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }
    numeric_level = level_map.get(level.lower())
    if numeric_level is None:
        raise ValueError(
            f"Invalid log level: {level!r}. "
            f"Must be one of: {', '.join(level_map.keys())}"
        )

    # Default format
    if format_string is None:
        if include_module_name:
            format_string = "{asctime:23} {levelname:>7s} | {name:35s} | {message}"
        else:
            format_string = "{asctime:23} {levelname:>7s} | {message}"

    formatter = logging.Formatter(format_string, style='{')

    # Get root logger for qhat package
    root_logger = logging.getLogger("qhat")
    root_logger.setLevel(numeric_level)

    # Clear any existing handlers (in case called multiple times)
    root_logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(numeric_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # File handler (if specified)
    if logfile:
        logfile_path = Path(logfile)
        logfile_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(logfile, mode='w')  # Overwrite mode
        file_handler.setLevel(numeric_level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

        root_logger.info(f"Logging to file: {logfile_path.absolute()}")

    # Don't propagate to root logger (keeps qhat logs separate from other packages)
    root_logger.propagate = False

    return root_logger


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger for a module.

    This is a convenience function. In most cases, prefer using
    `logging.getLogger(__name__)` directly in modules.

    Args:
        name: Module name (usually __name__)

    Returns:
        Logger instance

    Example:
        >>> from qhat.logging_utils import get_logger
        >>> logger = get_logger(__name__)
        >>> logger.info("Module loaded")
    """
    return logging.getLogger(name)


@contextmanager
def temporary_log_level(level: int, logger_name: str = "qhat"):
    """
    Context manager to temporarily change log level.

    Useful for making specific sections more or less verbose without
    affecting the rest of the application.

    Args:
        level: Numeric log level (e.g., logging.DEBUG, logging.INFO)
        logger_name: Name of logger to modify (default: "qhat" root logger)

    Yields:
        The logger with temporarily modified level

    Example:
        >>> import logging
        >>> from qhat.logging_utils import temporary_log_level
        >>>
        >>> # Temporarily enable debug logging
        >>> with temporary_log_level(logging.DEBUG, "qhat.analysis"):
        ...     # Verbose logging only in this block
        ...     result = some_computation()
    """
    logger = logging.getLogger(logger_name)
    old_level = logger.level
    try:
        logger.setLevel(level)
        yield logger
    finally:
        logger.setLevel(old_level)


@contextmanager
def log_execution_time(operation: str, logger: Optional[logging.Logger] = None):
    """
    Context manager to log execution time of an operation.

    Args:
        operation: Description of the operation being timed
        logger: Logger to use (uses qhat root logger if None)

    Example:
        >>> from qhat.logging_utils import log_execution_time
        >>> import logging
        >>> logger = logging.getLogger(__name__)
        >>>
        >>> with log_execution_time("matrix computation", logger):
        ...     result = expensive_computation()
        >>> # Logs: "matrix computation completed in 1.23s"
    """
    import time

    if logger is None:
        logger = logging.getLogger("qhat")

    logger.debug(f"Starting {operation}")
    start_time = time.time()
    try:
        yield
        elapsed = time.time() - start_time
        logger.info(f"{operation} completed in {elapsed:.2f}s")
    except Exception as e:
        elapsed = time.time() - start_time
        logger.error(f"{operation} failed after {elapsed:.2f}s: {e}")
        raise


def reconfigure_log_level(level: str):
    """
    Change the log level for all qhat loggers.

    This is useful for dynamically adjusting verbosity during execution.

    Args:
        level: New log level string ("debug", "verbose", "info", "warning", "error")

    Raises:
        ValueError: If log level is invalid
    """
    level_map = {
        "debug": logging.DEBUG,
        "verbose": VERBOSE,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }
    numeric_level = level_map.get(level.lower())
    if numeric_level is None:
        raise ValueError(
            f"Invalid log level: {level!r}. "
            f"Must be one of: {', '.join(level_map.keys())}"
        )

    # Update qhat root logger
    qhat_logger = logging.getLogger("qhat")
    qhat_logger.setLevel(numeric_level)

    # Update all handlers
    for handler in qhat_logger.handlers:
        handler.setLevel(numeric_level)
