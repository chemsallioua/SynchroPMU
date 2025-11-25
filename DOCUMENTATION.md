# Quick Documentation Guide

## Overview

This repository now includes comprehensive documentation generated with Doxygen. The documentation covers:

- ✅ **All public API functions** - Detailed descriptions with parameters and return values
- ✅ **All data structures** - Complete field-by-field documentation
- ✅ **Header files** - Full file-level documentation
- ✅ **Source implementation** - Algorithm overview and implementation details
- ✅ **Python API** - Complete docstrings for Python bindings
- ✅ **Usage examples** - Code examples in documentation comments

## Quick Links

### For API Users

1. **Main Documentation Entry Point**: Open `Doxygen/html/index.html` in your browser
2. **Function Reference**: Browse all functions at `Doxygen/html/globals_func.html`
3. **Data Structures**: View all structs at `Doxygen/html/annotated.html`
4. **Python API**: See docstrings in `python/pmu_estimator.py`

### For Developers

1. **Source Documentation**: `Doxygen/html/files.html` - All source files with detailed comments
2. **Function Stubs**: `Doxygen/html/func__stubs_8h.html` - Customization interface
3. **Algorithm Details**: See file header in `src/pmu_estimator.c`

## Key Documentation Features

### C API Documentation

#### Main Functions
- `pmu_init()` - Initialize PMU estimator with configuration
- `pmu_estimate()` - Estimate synchrophasor, frequency, and ROCOF
- `pmu_deinit()` - Clean up and free resources
- `pmu_dump_frame()` - Output results to stream

#### Key Structures
- `pmu_context` - Main PMU instance (supports multiple independent instances)
- `estimator_config` - Configuration parameters
- `pmu_frame` - Output frame with synchrophasor and ROCOF
- `phasor` - Synchrophasor representation (amplitude, phase, frequency)

### Python API Documentation

All Python classes and methods have detailed docstrings:
- `PMUEstimator` - Main estimator class
- `EstimatorConfig` - Configuration object
- `PmuFrame` - Output frame structure
- `Phasor` - Synchrophasor object

## Regenerating Documentation

After making changes to code comments:

```bash
doxygen Doxyfile
```

The updated HTML documentation will be generated in `Doxygen/html/`.

## Documentation Standards Used

- **Doxygen format** - Industry-standard documentation tool
- **Detailed parameter descriptions** - All input/output parameters documented
- **Return value documentation** - What each function returns
- **Usage examples** - Code snippets showing typical usage
- **Cross-references** - Links between related functions and structures
- **Platform notes** - Platform-specific information where applicable

## Finding Specific Information

### "How do I initialize the PMU?"
See `pmu_init()` documentation in `Doxygen/html/pmu__estimator_8h.html#a42a30dd06ee4f11afdd4af8850052507`

### "What configuration parameters are available?"
See `estimator_config` structure documentation in `Doxygen/html/structestimator__config.html`

### "How do I use the Python API?"
See the module docstring and class documentation in `python/pmu_estimator.py`

### "How does the algorithm work?"
See the algorithm overview in the file header of `src/pmu_estimator.c`

### "Can I customize mathematical functions?"
See `func_stubs.h` documentation for the function stubbing mechanism

## Coverage Summary

| Component | Documentation Status |
|-----------|---------------------|
| Public API Functions | ✅ Complete |
| Data Structures | ✅ Complete |
| Function Parameters | ✅ Complete |
| Return Values | ✅ Complete |
| Usage Examples | ✅ Included |
| Python API | ✅ Complete |
| Algorithm Overview | ✅ Complete |
| Configuration Options | ✅ Complete |
| Build Instructions | ✅ In README |

## Next Steps

1. **Browse the documentation**: Start at `Doxygen/html/index.html`
2. **Read the main README**: See `README.md` for build and usage instructions
3. **Try the examples**: Check `examples/` directory for working code
4. **Explore the Python API**: See `python/` for Python bindings

For questions or issues, refer to the detailed function documentation or the source code comments.
