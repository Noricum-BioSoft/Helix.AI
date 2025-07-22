# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased] - 2024-01-XX

### 🎨 **UI/UX Improvements**
- **Enhanced Layout**: Implemented professional 75/25 layout with optimized spacing
- **Tips Section**: Added helpful tips section positioned above command input for better UX
- **Clean Header**: Simplified header to just "Helix.AI" for cleaner appearance
- **Improved Spacing**: Balanced margins and padding (10-15% whitespace) for professional look
- **File Upload**: Removed auto-population on drag-and-drop for cleaner user experience
- **Placeholder Text**: Updated command input placeholder to "visualize the phylogenetic tree"

### 🔬 **New Features**
- **Clustering Analysis**: Added hierarchical clustering with representative sequence selection
- **Average Distance Calculation**: Implemented sophisticated distance metrics for cluster analysis
- **Plasmid Visualization**: Enhanced with circular/linear/both view options using SeqViz
- **Command Routing**: Improved natural language command routing to appropriate tools
- **Workflow Context**: Enhanced data passing between different analysis steps

### 🐛 **Bug Fixes**
- **Command Routing**: Fixed plasmid visualization commands being incorrectly routed to alignment
- **Import Errors**: Resolved missing component imports in frontend
- **JSON Serialization**: Fixed numpy type serialization issues in backend responses
- **Session Management**: Improved session context handling and persistence

### 📚 **Documentation**
- **Updated README**: Comprehensive documentation of all features and improvements
- **Code Comments**: Enhanced code documentation and inline comments
- **API Documentation**: Updated endpoint documentation and examples
- **User Guides**: Improved natural language command examples and workflows

### 🧪 **Testing**
- **Test Coverage**: Enhanced test suite for new clustering and plasmid features
- **Integration Tests**: Added end-to-end workflow testing
- **Performance Tests**: Improved backend performance monitoring

## [0.2.0] - 2024-01-XX

### 🌳 **Phylogenetic Tree Enhancements**
- **ETE3 Integration**: Added high-quality phylogenetic tree visualization with SVG rendering
- **Tree Statistics**: Comprehensive tree metrics and relationship analysis
- **Interactive Features**: Zoom, pan, and interactive tree exploration
- **Newick Support**: Full support for standard Newick format trees

### 🔬 **DNA Synthesis Vendor Research**
- **Vendor Comparison**: Compare 7 major DNA synthesis vendors
- **Pricing Analysis**: Real-time pricing ranges and service comparisons
- **Service Matching**: Find vendors based on sequence length and quantity
- **AI Recommendations**: Intelligent vendor selection advice

### 🤖 **Natural Language Commands**
- **Enhanced Parsing**: Improved natural language command understanding
- **Multi-step Workflows**: Chain commands for complex bioinformatics pipelines
- **Context Awareness**: Maintains session history for continuous workflows
- **Smart Routing**: Intelligent routing to appropriate bioinformatics tools

### 📊 **Visualization Improvements**
- **SeqViz Integration**: Professional plasmid and vector visualization
- **Plotly Charts**: Dynamic charts for sequence analysis and mutations
- **Real-time Updates**: Live visualizations that update as you work
- **Export Capabilities**: Save plots and results in multiple formats

### 🧰 **Tool Suite Enhancements**
- **Sequence Alignment**: Multiple algorithms (ClustalW, Muscle, MAFFT)
- **Mutation Analysis**: Generate and analyze sequence variants
- **Data Science Tools**: Statistical analysis and feature engineering
- **Variant Selection**: Smart selection based on diversity and custom criteria

### 🔄 **Session Management**
- **Persistent Sessions**: Track workflow across multiple commands
- **History Tracking**: Complete audit trail of all operations
- **Context Preservation**: Maintain state between commands
- **Workflow Context**: Pass data between different analysis steps

### 🎯 **User Experience**
- **Drag-and-Drop**: Upload FASTA/CSV files directly
- **Command Mode Toggle**: Switch between natural language and structured commands
- **Real-time Feedback**: Immediate response and progress indicators
- **Responsive Design**: Works on desktop and mobile devices

## [0.1.0] - 2024-01-XX

### 🚀 **Initial Release**
- **FastAPI Backend**: Modern web framework with MCP integration
- **React Frontend**: TypeScript-based user interface
- **Natural Language Processing**: Basic command parsing and execution
- **Session Management**: Basic session tracking and history
- **File Upload**: Drag-and-drop file upload functionality
- **Basic Tools**: Initial set of bioinformatics tools
- **Documentation**: Basic setup and usage documentation

### 🔧 **Core Features**
- **Sequence Alignment**: Basic sequence alignment capabilities
- **Mutation Generation**: Simple sequence mutation tools
- **Data Analysis**: Basic statistical analysis features
- **Visualization**: Initial plotting and chart capabilities
- **API Endpoints**: RESTful API for tool execution
- **Health Monitoring**: Basic health check endpoints

### 📚 **Documentation**
- **Setup Guide**: Installation and configuration instructions
- **API Documentation**: Endpoint documentation and examples
- **User Guide**: Basic usage instructions
- **Development Guide**: Contributing guidelines

---

## Version History

- **0.1.0**: Initial release with basic bioinformatics tools
- **0.2.0**: Enhanced phylogenetic analysis and vendor research
- **Unreleased**: UI improvements, clustering analysis, and comprehensive documentation

## Contributing

When contributing to this project, please update this changelog with a brief description of your changes under the appropriate version section. 