# TurboRVB Documentation

This directory contains the automatically generated documentation for the TurboRVB project using Doxygen.

## Overview

TurboRVB is a Quantum Monte Carlo package for electronic structure calculations. This documentation provides comprehensive API reference, algorithm descriptions, and usage examples for all Fortran subroutines and functions in the codebase.

## Documentation Structure

- **HTML Documentation**: `html/index.html` - Interactive web-based documentation
- **LaTeX Documentation**: `latex/refman.pdf` - PDF documentation for printing
- **Source Code**: Cross-referenced source code with syntax highlighting

## Generating Documentation

### Prerequisites

1. **Install Doxygen**:
   ```bash
   # macOS
   brew install doxygen
   
   # Ubuntu/Debian
   sudo apt-get install doxygen
   
   # CentOS/RHEL
   sudo yum install doxygen
   
   # Windows (using Chocolatey)
   choco install doxygen
   ```

2. **Optional Dependencies**:
   ```bash
   # For PDF generation (LaTeX)
   # macOS
   brew install --cask mactex
   
   # Ubuntu/Debian
   sudo apt-get install texlive-full
   
   # For Graphviz diagrams (if enabled)
   # macOS
   brew install graphviz
   
   # Ubuntu/Debian
   sudo apt-get install graphviz
   ```

### Quick Start

1. **Navigate to project root**:
   ```bash
   cd /path/to/turborvb
   ```

2. **Generate documentation**:
   ```bash
   doxygen Doxyfile
   ```

3. **View documentation**:
   ```bash
   # Open HTML documentation
   open docs/html/index.html
   
   # Or use a web browser
   firefox docs/html/index.html
   ```

### Advanced Usage

#### Custom Configuration

The documentation is generated using the `Doxyfile` configuration. You can modify settings by:

1. **Editing Doxyfile directly**:
   ```bash
   # Edit the configuration
   nano Doxyfile
   
   # Regenerate documentation
   doxygen Doxyfile
   ```

2. **Creating a custom configuration**:
   ```bash
   # Generate a template
   doxygen -g my_doxyfile
   
   # Edit the template
   nano my_doxyfile
   
   # Use custom configuration
   doxygen my_doxyfile
   ```

#### Common Configuration Modifications

- **Change output directory**: Modify `OUTPUT_DIRECTORY` in Doxyfile
- **Include/exclude files**: Modify `INPUT` and `EXCLUDE` patterns
- **Enable/disable features**: Modify `GENERATE_HTML`, `GENERATE_LATEX`, etc.
- **Customize appearance**: Modify `HTML_EXTRA_STYLESHEET` for custom CSS

### Continuous Integration

For automated documentation generation in CI/CD pipelines:

```yaml
# Example GitHub Actions workflow
name: Generate Documentation
on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  docs:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Install Doxygen
      run: sudo apt-get install doxygen
    - name: Generate Documentation
      run: doxygen Doxyfile
    - name: Deploy to GitHub Pages
      uses: peaceiris/actions-gh-pages@v3
      with:
        github_token: ${{ secrets.GITHUB_TOKEN }}
        publish_dir: ./docs/html
```

## Documentation Features

### HTML Documentation
- **Interactive Navigation**: Tree view and search functionality
- **Cross-references**: Links between functions, modules, and files
- **Source Code**: Syntax-highlighted source code with line numbers
- **Search Engine**: Full-text search across all documentation
- **Responsive Design**: Works on desktop and mobile devices

### LaTeX/PDF Documentation
- **Print-friendly**: Optimized for printing and offline reading
- **Hyperlinks**: Clickable references and cross-references
- **Table of Contents**: Comprehensive index and navigation
- **Mathematical Formulas**: Proper rendering of LaTeX equations

### Source Code Integration
- **Inline Documentation**: Function bodies included in documentation
- **Parameter Documentation**: Detailed descriptions of all parameters
- **Return Value Documentation**: Clear explanation of function outputs
- **Usage Examples**: Code examples and usage patterns

## Troubleshooting

### Common Issues

1. **Doxygen not found**:
   ```bash
   # Check installation
   which doxygen
   
   # Reinstall if needed
   brew reinstall doxygen  # macOS
   sudo apt-get install --reinstall doxygen  # Ubuntu
   ```

2. **LaTeX compilation errors**:
   ```bash
   # Install missing LaTeX packages
   sudo apt-get install texlive-latex-extra
   
   # Or disable LaTeX generation
   # Set GENERATE_LATEX = NO in Doxyfile
   ```

3. **Memory issues with large projects**:
   ```bash
   # Increase memory limit
   export DOXYGEN_MEMORY_LIMIT=2048
   doxygen Doxyfile
   ```

4. **Missing source files**:
   ```bash
   # Check INPUT configuration in Doxyfile
   # Ensure all source directories are included
   ```

### Performance Optimization

For large codebases, consider:

1. **Parallel processing**:
   ```bash
   # Set NUM_PROC_THREADS in Doxyfile
   NUM_PROC_THREADS = 4
   ```

2. **Selective documentation**:
   ```bash
   # Use EXCLUDE patterns to skip certain files
   EXCLUDE = */test/* */examples/*
   ```

3. **Incremental generation**:
   ```bash
   # Only regenerate changed files
   doxygen Doxyfile
   ```

## Customization

### Adding Custom Documentation

1. **Module Documentation**:
   ```fortran
   !> @brief Brief description of the module
   !> @details Detailed description of functionality
   module my_module
   ```

2. **Function Documentation**:
   ```fortran
   !> @brief Brief description of the function
   !> @param[in] input_param Description of input parameter
   !> @param[out] output_param Description of output parameter
   !> @return Description of return value
   function my_function(input_param, output_param)
   ```

3. **Algorithm Documentation**:
   ```fortran
   !> @details
   !> This subroutine implements the following algorithm:
   !> 1. Step one description
   !> 2. Step two description
   !> 3. Step three description
   !>
   !> Mathematical formulation:
   !> \[ E = \sum_{i} \langle \psi_i | H | \psi_i \rangle \]
   subroutine my_algorithm()
   ```

### Styling and Branding

1. **Custom CSS**:
   ```bash
   # Create custom stylesheet
   echo "body { font-family: 'Arial', sans-serif; }" > custom.css
   
   # Add to Doxyfile
   HTML_EXTRA_STYLESHEET = custom.css
   ```

2. **Custom Header/Footer**:
   ```bash
   # Create custom header
   doxygen -w html header.html footer.html stylesheet.css
   
   # Edit and use in Doxyfile
   HTML_HEADER = header.html
   HTML_FOOTER = footer.html
   ```

## Contributing to Documentation

### Documentation Standards

1. **Use English**: All documentation should be in English
2. **Be Consistent**: Follow established patterns and conventions
3. **Include Examples**: Provide usage examples where appropriate
4. **Document Parameters**: All parameters should have clear descriptions
5. **Explain Algorithms**: Mathematical algorithms should be well-documented

### Adding New Documentation

1. **Follow Doxygen Standards**: Use proper Doxygen tags and syntax
2. **Test Generation**: Always test documentation generation after changes
3. **Update This README**: Keep this file updated with new features or changes

## Support

For issues with documentation generation:

1. **Check Doxygen Version**: Ensure you're using a compatible version
2. **Review Configuration**: Verify Doxyfile settings
3. **Check Dependencies**: Ensure all required tools are installed
4. **Consult Doxygen Manual**: [https://www.doxygen.nl/manual/](https://www.doxygen.nl/manual/)

## License

This documentation is generated from the TurboRVB source code and follows the same license as the project (GNU General Public License v3.0).

---

**Note**: This documentation is automatically generated from source code comments. For the most up-to-date information, always refer to the source code and generated documentation. 