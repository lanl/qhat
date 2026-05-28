"""
Simple test script to verify the new logging architecture works correctly.

This script demonstrates:
1. Configuring logging at startup
2. Using module-level loggers
3. Different log levels (debug, verbose, info, warning, error)
4. Temporary log level changes
5. Execution time logging
"""

import logging
from qhat.logging_utils import configure_logging, temporary_log_level, log_execution_time
import time

# Configure logging
configure_logging(level="verbose", logfile="test_logging.log")

# Get module logger
logger = logging.getLogger(__name__)

def test_basic_logging():
    """Test basic logging at different levels."""
    logger.debug("This is a DEBUG message")
    logger.verbose("This is a VERBOSE message")
    logger.info("This is an INFO message")
    logger.warning("This is a WARNING message")
    logger.error("This is an ERROR message")

def test_temporary_level():
    """Test temporary log level changes."""
    logger.info("Normal verbosity level")

    with temporary_log_level(logging.DEBUG):
        logger.debug("This DEBUG message is visible inside the context manager")

    logger.debug("This DEBUG message is NOT visible outside")
    logger.info("Back to normal verbosity")

def test_execution_time():
    """Test execution time logging."""
    with log_execution_time("test operation", logger):
        time.sleep(0.1)
        logger.info("Doing some work...")

def test_module_specific_logging():
    """Test that different modules can have different loggers."""
    analysis_logger = logging.getLogger("qhat.analysis")
    common_logger = logging.getLogger("qhat.common")

    analysis_logger.info("Message from analysis module")
    common_logger.info("Message from common module")

if __name__ == "__main__":
    logger.info("=" * 80)
    logger.info("Testing QHAT Logging Architecture")
    logger.info("=" * 80)

    logger.info("\n1. Testing basic logging levels:")
    test_basic_logging()

    logger.info("\n2. Testing temporary log level:")
    test_temporary_level()

    logger.info("\n3. Testing execution time logging:")
    test_execution_time()

    logger.info("\n4. Testing module-specific logging:")
    test_module_specific_logging()

    logger.info("\n" + "=" * 80)
    logger.info("Logging test complete!")
    logger.info("=" * 80)
