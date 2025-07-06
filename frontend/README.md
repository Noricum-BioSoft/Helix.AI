# DataBloom.AI Frontend

This frontend application provides a user-friendly interface for interacting with the Bioinformatics MCP Server.

## Features

- **Intelligent Command Parsing**: Automatically routes user commands to appropriate MCP endpoints
- **Real-time Server Status**: Shows connection status to the backend MCP server
- **Available Tools Display**: Lists all available MCP tools with descriptions and parameters
- **Rich Output Rendering**: Displays results in various formats (tables, plots, text)
- **Command History**: Maintains a history of all executed commands with timestamps

## Architecture

### Components

- **App.tsx**: Main application component with command interface
- **services/mcpApi.ts**: API service for communicating with MCP server
- **utils/commandParser.ts**: Intelligent command parsing and routing

### Command Routing

The application uses intelligent command parsing to route user input to the appropriate MCP endpoints:

1. **Sequence Alignment**: Commands containing "align", "alignment", "clustal", etc.
2. **Mutation Analysis**: Commands containing "mutate", "variant", "mutation", etc.
3. **Data Analysis**: Commands containing "analyze", "analysis", "phylogeny", etc.
4. **Visualization**: Commands containing "visualize", "plot", "chart", etc.
5. **General Commands**: Fallback to the legacy `/execute` endpoint

## Installation

1. Install dependencies:
```bash
cd frontend
npm install
```

2. Start the development server:
```bash
npm run dev
```

3. Ensure the backend MCP server is running:
```bash
cd ../backend
python main_with_mcp.py
```

## Usage

### Example Commands

#### Sequence Alignment
```
align sequences ACTGTTGAC ACTGCATCC
align with clustal algorithm
multiple sequence alignment
```

#### Mutation Analysis
```
mutate sequence ACTGTTGAC
generate 10 variants of ACTGTTGAC
create mutations with 0.2 mutation rate
```

#### Data Analysis
```
analyze sequence data
analyze for phylogeny
sequence composition analysis
```

#### Visualization
```
visualize alignment
plot alignment in PNG format
create alignment chart
```

### Interface Features

1. **Server Status Indicator**: Shows if the backend is healthy
2. **Available Tools Sidebar**: Lists all MCP tools with parameters
3. **Command History**: Shows all executed commands with results
4. **Rich Output Display**: Renders results in appropriate formats

## API Integration

The frontend communicates with the backend through these endpoints:

- `POST /execute` - General command execution (legacy)
- `POST /mcp/sequence-alignment` - Sequence alignment
- `POST /mcp/mutate-sequence` - Sequence mutation
- `POST /mcp/analyze-sequence-data` - Data analysis
- `POST /mcp/visualize-alignment` - Alignment visualization
- `GET /mcp/tools` - List available tools
- `GET /health` - Server health check

## Development

### Adding New Command Types

1. Update `CommandParser` in `utils/commandParser.ts`:
   - Add new keywords to the appropriate keyword arrays
   - Create a new parsing method
   - Update the main `parseCommand` method

2. Update the API service in `services/mcpApi.ts`:
   - Add new interface for the request
   - Add new method to call the endpoint

3. Update the App component in `App.tsx`:
   - Add new case to the switch statement in `handleSubmit`
   - Update the `ParsedCommand` interface

### Styling

The application uses Bootstrap 5 for styling. Custom styles can be added to:
- `src/App.css` for component-specific styles
- `src/index.css` for global styles

### Error Handling

The application includes comprehensive error handling:
- Network errors are displayed to the user
- Invalid commands show helpful error messages
- Server errors are clearly indicated
- Loading states prevent multiple submissions

## Troubleshooting

### Common Issues

1. **Server Connection Failed**:
   - Ensure the backend is running on port 8001
   - Check CORS settings in the backend
   - Verify the API_BASE_URL in mcpApi.ts

2. **Commands Not Working**:
   - Check the command parser keywords
   - Verify the MCP server has the required tools
   - Check browser console for errors

3. **TypeScript Errors**:
   - Run `npm run build` to check for type errors
   - Ensure all interfaces are properly defined
   - Check import/export statements

### Debugging

1. Open browser developer tools
2. Check the Network tab for API calls
3. Check the Console tab for errors
4. Use the React Developer Tools extension

## Building for Production

```bash
npm run build
```

The built files will be in the `dist/` directory and can be served by any static file server.

## Contributing

When adding new features:

1. Follow the existing code structure
2. Add proper TypeScript types
3. Include error handling
4. Update this README with new features
5. Test with various command inputs 