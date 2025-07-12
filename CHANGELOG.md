# Changelog

All notable changes to the DataBloom.AI project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **🌳 Interactive Phylogenetic Tree Visualization**
  - D3.js-based interactive tree rendering with zoom and pan
  - Newick format parsing and display
  - Comprehensive tree statistics and metrics
  - Real-time tree updates during analysis
  - Support for complex tree structures with proper error handling

- **🔬 DNA Synthesis Vendor Research**
  - Complete vendor comparison system with 7 major vendors
  - Real-time pricing analysis and service comparisons
  - AI-powered vendor recommendations based on sequence requirements
  - Interactive vendor cards with detailed specifications
  - Service matching based on sequence length and quantity

- **Natural Language Command Processing**
  - Intelligent command parser that understands natural language bioinformatics commands
  - Support for complex multi-step workflows
  - Context-aware command execution with session management
  - Example: "visualize the variants in a phylogenetic tree"

- **Session Management System**
  - Persistent session tracking across multiple commands
  - Complete audit trail of all operations
  - Session-based history management
  - Automatic session creation and cleanup

- **Enhanced Frontend Interface**
  - Command mode toggle (Natural Language vs Structured Commands)
  - Drag-and-drop file upload support for FASTA and CSV files
  - Real-time feedback and progress indicators
  - Responsive design for desktop and mobile

- **Comprehensive Tool Suite**
  - **Mutation Tool**: Generate and analyze sequence variants with statistics
  - **Alignment Tool**: Multiple sequence alignment with ClustalW/Muscle/MAFFT
  - **Data Science Tool**: Statistical analysis, correlation matrices, feature engineering
  - **Variant Selection Tool**: Smart selection based on diversity, length, or custom criteria
  - **Phylogenetic Tree Tool**: Tree construction and analysis with UPGMA method
  - **DNA Vendor Research Tool**: Vendor comparison and recommendations

- **Visualization Enhancements**
  - Plotly integration for dynamic charts
  - Real-time plot updates
  - Export capabilities for plots and results
  - Interactive visualizations for sequence analysis
  - D3.js phylogenetic tree visualization

### Changed
- **Backend Architecture**
  - Migrated from basic MCP server to enhanced FastAPI with MCP integration
  - Improved error handling and validation
  - Better dependency management with conda environment
  - Enhanced API endpoints with session support
  - Updated command router for better tool integration

- **Frontend Architecture**
  - Updated to use axios for API calls
  - Improved state management for session handling
  - Enhanced UI with Bootstrap components
  - Better error handling and user feedback
  - Added D3.js integration for phylogenetic trees
  - Enhanced vendor research display components

### Fixed
- **Import and Dependency Issues**
  - Resolved Python path issues for tool imports
  - Fixed conda environment compatibility
  - Updated dependency versions for better stability
  - Fixed D3.js TypeScript integration issues

- **Session Data Protection**
  - Added comprehensive `.gitignore` patterns to prevent session data from being committed
  - Removed tracked session files from git
  - Added patterns for user data protection

- **Phylogenetic Tree Rendering**
  - Fixed NaN coordinate issues in D3.js tree rendering
  - Improved Newick parser for complex tree structures
  - Enhanced error handling for malformed tree data
  - Fixed tree layout and positioning issues

- **Vendor Research Display**
  - Fixed nested data structure handling in frontend
  - Improved vendor card layout and information display
  - Enhanced recommendations display
  - Fixed API response parsing for vendor data

## [1.0.0] - 2024-01-XX

### Added
- **Initial MCP Server Implementation**
  - Basic Model Context Protocol server
  - Standard bioinformatics tool integration
  - JSON-RPC communication protocol

- **React Frontend**
  - Bootstrap-based responsive UI
  - Command input interface
  - Result display with basic formatting

- **Core Bioinformatics Tools**
  - Sequence alignment capabilities
  - Basic mutation generation
  - Sequence validation and analysis

- **FastAPI Backend**
  - RESTful API endpoints
  - CORS middleware for frontend integration
  - Health check and status endpoints

### Changed
- **Project Structure**
  - Organized into frontend/backend/tools structure
  - Modular tool system for easy extension
  - Clear separation of concerns

### Fixed
- **Initial Setup Issues**
  - Resolved Python environment setup
  - Fixed Node.js dependency installation
  - Corrected import paths and module resolution

---

## Version History

### [0.9.0] - Development Phase
- Initial project setup
- Basic MCP server implementation
- React frontend foundation
- Core bioinformatics tool development

### [0.8.0] - Alpha Release
- Enhanced MCP server with validation
- Improved frontend interface
- Better error handling
- Documentation updates

### [0.7.0] - Beta Release
- Natural language command processing
- Session management system
- Comprehensive tool suite
- Advanced visualizations

---

## Migration Guide

### From Version 0.8.0 to 1.0.0

#### Backend Changes
1. **Environment Setup**
   ```bash
   # Update to Python 3.10+
   conda create -n databloom python=3.10
   conda activate databloom
   pip install -r requirements.txt
   ```

2. **New Dependencies**
   ```bash
   # Install additional packages
   conda install plotly pandas seaborn
   pip install langchain-deepseek
   ```

3. **Configuration Updates**
   ```bash
   # Update environment variables
   export SESSION_TIMEOUT=3600
   export MAX_SESSION_SIZE=1000
   ```

#### Frontend Changes
1. **New Dependencies**
   ```bash
   npm install axios @types/axios
   ```

2. **API Service Updates**
   - Updated to use new session-aware endpoints
   - Added natural language command support
   - Enhanced error handling

#### Data Migration
- No data migration required
- Session data is automatically managed
- Old command format still supported

---

## Breaking Changes

### Version 1.0.0
- **API Endpoint Changes**
  - New session management endpoints
  - Updated natural language command endpoints
  - Enhanced response format for better frontend integration

- **Environment Requirements**
  - Python 3.10+ required (up from 3.9)
  - Conda environment recommended
  - Additional bioinformatics tools required

- **Frontend Requirements**
  - Node.js 16+ required
  - Additional npm packages required
  - Updated API service integration

---

## Known Issues

### Version 1.0.0
- **Session Management**
  - Large session files may impact performance
  - Session cleanup not yet automated
  - Workaround: Manual session cleanup recommended

- **Natural Language Processing**
  - Limited to English commands
  - Complex nested commands may not parse correctly
  - Workaround: Use structured commands for complex operations

- **Visualization**
  - Plotly charts may not render in some browsers
  - Large datasets may cause performance issues
  - Workaround: Use smaller datasets or export to static images

---

## Future Roadmap

### Version 1.1.0 (Planned)
- **Enhanced Natural Language Processing**
  - Multi-language support
  - More complex command parsing
  - Better context understanding

- **Advanced Analytics**
  - Machine learning integration
  - Predictive modeling
  - Advanced statistical analysis

- **Performance Improvements**
  - Caching system for repeated operations
  - Background processing for long-running tasks
  - Optimized memory usage

### Version 1.2.0 (Planned)
- **Collaboration Features**
  - Multi-user sessions
  - Shared workflows
  - Real-time collaboration

- **Advanced Visualization**
  - 3D molecular visualization
  - Interactive phylogenetic trees
  - Custom plot templates

- **Integration Features**
  - External database connections
  - API integrations with other bioinformatics tools
  - Export to various formats

---

## Contributing

When contributing to this project, please:

1. **Update the CHANGELOG** for any user-facing changes
2. **Follow the existing format** for entries
3. **Include migration notes** for breaking changes
4. **Document new features** with examples
5. **Update version numbers** according to semantic versioning

### Changelog Entry Format

```markdown
### Added
- New feature description
- Another new feature

### Changed
- Changed feature description

### Fixed
- Bug fix description

### Removed
- Removed feature description
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 