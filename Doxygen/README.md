# SynchroPMU Documentation

This directory contains the comprehensive API documentation for the SynchroPMU library, generated using Doxygen.

## Contents

- **html/** - HTML documentation (browse starting from `index.html`)

## Viewing the Documentation

### Online

If this repository is hosted on GitHub, you can view the documentation by:
1. Opening the `html/index.html` file in your web browser
2. Or deploying the `html/` directory to GitHub Pages

### Local

To view the documentation locally:

```bash
# Navigate to the Doxygen/html directory
cd Doxygen/html

# Open index.html in your default browser
# On Linux:
xdg-open index.html

# On macOS:
open index.html

# On Windows:
start index.html
```

## Regenerating Documentation

If you make changes to the source code documentation comments, you can regenerate the documentation:

```bash
# From the repository root directory
doxygen Doxyfile
```

This will update the HTML documentation in the `Doxygen/html/` directory.

## Documentation Structure

The documentation includes:

### API Reference
- **Data Types and Structures**: All structs, typedefs, and enumerations
- **Functions**: Detailed descriptions of all public API functions
- **Macros**: Preprocessor macros and their usage
- **Files**: Source file documentation with dependencies

### Key Sections
- **Main Page**: Overview of the library (from README.md)
- **Data Structures**: Complete list of structures with field descriptions
- **File List**: All documented source files
- **Function Index**: Alphabetical list of all functions

### Main Components Documented

1. **PMU Estimator Core** (`pmu_estimator.h`, `pmu_estimator.c`)
   - PMU initialization and deinitialization
   - Synchrophasor estimation
   - ROCOF calculation
   - Configuration management

2. **Function Stubs** (`func_stubs.h`)
   - Type definitions
   - Mathematical function stubs
   - Customization interface

3. **Data Structures**
   - `pmu_context`: Main PMU instance context
   - `estimator_config`: Configuration parameters
   - `phasor`: Synchrophasor representation
   - `pmu_frame`: Output frame structure

## Requirements for Regeneration

To regenerate the documentation yourself, you need:
- Doxygen (version 1.9.8 or higher)
- Graphviz (for generating diagrams)

Install on Ubuntu/Debian:
```bash
sudo apt-get install doxygen graphviz
```

Install on macOS:
```bash
brew install doxygen graphviz
```

## Configuration

The Doxygen configuration is stored in `Doxyfile` in the repository root. Key settings:

- **PROJECT_NAME**: SynchroPMU
- **INPUT**: `src/` and `README.md`
- **OUTPUT_DIRECTORY**: `Doxygen`
- **RECURSIVE**: Yes
- **EXTRACT_ALL**: Yes
- **HAVE_DOT**: Yes (for graphs)

## Contributing

When adding new functions or modifying existing ones:

1. Add Doxygen-style comments to your code:
   ```c
   /**
    * @brief Brief description
    * 
    * Detailed description here.
    * 
    * @param[in] param1 Description of param1
    * @param[out] param2 Description of param2
    * @return Return value description
    */
   ```

2. Regenerate the documentation:
   ```bash
   doxygen Doxyfile
   ```

3. Review the generated documentation to ensure accuracy

## License

The documentation follows the same license as the SynchroPMU library.

Copyright (c) 2023. All Rights Reserved.
Confidential and Proprietary - University of Bologna.
