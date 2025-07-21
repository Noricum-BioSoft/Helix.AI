"""
Directed Evolution Command Handler

This module provides command handling for directed evolution workflows,
integrating with the platform's command routing system.
"""

import json
import sys
import os
from typing import Dict, List, Any, Optional
from datetime import datetime

# Add the workflow engine directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'workflow-engine'))

from directed_evolution_workflow import DirectedEvolutionWorkflow, create_workflow_config


class DirectedEvolutionHandler:
    """Handler for directed evolution commands."""
    
    def __init__(self):
        self.workflow = DirectedEvolutionWorkflow()
        self.active_workflows = {}
    
    def handle_start_evolution(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle start evolution command."""
        try:
            # Extract parameters
            target_property = params.get('target_property', 'thermal_stability')
            library_size = params.get('library_size', 50)
            strategy = params.get('strategy', 'rational')
            assay_type = params.get('assay_type', 'thermal_stability')
            protein_name = params.get('name', 'Directed_Evolution_Protein')
            protein_description = params.get('description', 'Protein for directed evolution optimization')
            sequence_length = params.get('length', 250)
            
            # Create workflow configuration
            config = {
                'target_property': target_property,
                'library_size': library_size,
                'strategy': strategy,
                'assay_type': assay_type,
                'name': protein_name,
                'description': protein_description,
                'length': sequence_length,
                'import_sequence': False
            }
            
            # Run workflow
            results = self.workflow.run_complete_workflow(config)
            
            # Store workflow results
            workflow_id = results['workflow_info']['workflow_id']
            self.active_workflows[workflow_id] = results
            
            return {
                'status': 'success',
                'message': 'Directed evolution workflow completed successfully',
                'workflow_id': workflow_id,
                'results': {
                    'best_mutant': results['results']['best_mutant'].name,
                    'improvement': results['results']['report']['assay_results']['improvement'],
                    'total_mutants': results['results']['library'].library_size,
                    'target_property': target_property,
                    'assay_type': assay_type
                },
                'workflow_info': results['workflow_info']
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Failed to start evolution workflow: {str(e)}',
                'error': str(e)
            }
    
    def handle_import_sequence_evolution(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle evolution with imported sequence."""
        try:
            # Extract parameters
            sequence = params.get('sequence')
            if not sequence:
                return {
                    'status': 'error',
                    'message': 'No sequence provided for import'
                }
            
            target_property = params.get('target_property', 'thermal_stability')
            library_size = params.get('library_size', 50)
            strategy = params.get('strategy', 'rational')
            assay_type = params.get('assay_type', 'thermal_stability')
            protein_name = params.get('name', 'Imported_Protein')
            protein_description = params.get('description', 'Imported protein for evolution')
            
            # Create workflow configuration
            config = {
                'target_property': target_property,
                'library_size': library_size,
                'strategy': strategy,
                'assay_type': assay_type,
                'name': protein_name,
                'description': protein_description,
                'import_sequence': True,
                'sequence': sequence
            }
            
            # Run workflow
            results = self.workflow.run_complete_workflow(config)
            
            # Store workflow results
            workflow_id = results['workflow_info']['workflow_id']
            self.active_workflows[workflow_id] = results
            
            return {
                'status': 'success',
                'message': 'Directed evolution with imported sequence completed',
                'workflow_id': workflow_id,
                'results': {
                    'best_mutant': results['results']['best_mutant'].name,
                    'improvement': results['results']['report']['assay_results']['improvement'],
                    'total_mutants': results['results']['library'].library_size,
                    'target_property': target_property,
                    'assay_type': assay_type
                },
                'workflow_info': results['workflow_info']
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Failed to run evolution with imported sequence: {str(e)}',
                'error': str(e)
            }
    
    def handle_multi_cycle_evolution(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle multi-cycle evolution."""
        try:
            # Extract parameters
            num_cycles = params.get('num_cycles', 3)
            target_property = params.get('target_property', 'thermal_stability')
            library_size = params.get('library_size', 30)
            strategy = params.get('strategy', 'rational')
            assay_type = params.get('assay_type', 'thermal_stability')
            
            cycles = []
            current_sequence = None
            
            for cycle_num in range(1, num_cycles + 1):
                print(f"\n🔄 Running Evolution Cycle {cycle_num}/{num_cycles}")
                
                # Create workflow configuration for this cycle
                config = {
                    'target_property': target_property,
                    'library_size': library_size,
                    'strategy': strategy,
                    'assay_type': assay_type,
                    'name': f'Cycle_{cycle_num}_Protein',
                    'description': f'Protein for evolution cycle {cycle_num}',
                    'length': 250,
                    'import_sequence': False
                }
                
                # If we have a previous best mutant, use it as starting point
                if current_sequence:
                    config['import_sequence'] = True
                    config['sequence'] = current_sequence
                    config['name'] = f'Cycle_{cycle_num}_Mutant'
                
                # Run workflow
                workflow = DirectedEvolutionWorkflow()
                results = workflow.run_complete_workflow(config)
                
                # Store cycle results
                cycle_info = {
                    'cycle_number': cycle_num,
                    'workflow_id': results['workflow_info']['workflow_id'],
                    'best_mutant': results['results']['best_mutant'],
                    'improvement': results['results']['report']['assay_results']['improvement'],
                    'total_mutants': results['results']['library'].library_size,
                    'results': results
                }
                cycles.append(cycle_info)
                
                # Update current sequence for next cycle
                current_sequence = results['results']['best_mutant'].sequence
                
                # Reduce library size for subsequent cycles
                library_size = max(15, library_size - 5)
            
            # Calculate overall improvement
            initial_value = cycles[0]['results']['results']['sequence'].properties.get('thermal_stability', 0)
            final_value = cycles[-1]['results']['results']['best_mutant'].properties['thermal_stability']
            total_improvement = final_value - initial_value
            
            return {
                'status': 'success',
                'message': f'Multi-cycle evolution completed successfully',
                'cycles': cycles,
                'summary': {
                    'total_cycles': num_cycles,
                    'initial_value': initial_value,
                    'final_value': final_value,
                    'total_improvement': total_improvement,
                    'improvement_per_cycle': total_improvement / num_cycles,
                    'best_final_mutant': cycles[-1]['best_mutant'].name
                }
            }
            
        except Exception as e:
            return {
                'status': 'error',
                'message': f'Failed to run multi-cycle evolution: {str(e)}',
                'error': str(e)
            }
    
    def handle_get_workflow_status(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle get workflow status command."""
        workflow_id = params.get('workflow_id')
        
        if not workflow_id:
            return {
                'status': 'error',
                'message': 'No workflow ID provided'
            }
        
        if workflow_id not in self.active_workflows:
            return {
                'status': 'error',
                'message': f'Workflow {workflow_id} not found'
            }
        
        workflow_data = self.active_workflows[workflow_id]
        
        return {
            'status': 'success',
            'workflow_id': workflow_id,
            'workflow_info': workflow_data['workflow_info'],
            'steps': workflow_data['steps'],
            'results': {
                'best_mutant': workflow_data['results']['best_mutant'].name,
                'improvement': workflow_data['results']['report']['assay_results']['improvement'],
                'total_mutants': workflow_data['results']['library'].library_size
            }
        }
    
    def handle_list_workflows(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle list workflows command."""
        workflows = []
        
        for workflow_id, workflow_data in self.active_workflows.items():
            workflows.append({
                'workflow_id': workflow_id,
                'status': workflow_data['workflow_info']['status'],
                'start_time': workflow_data['workflow_info']['start_time'],
                'total_steps': workflow_data['workflow_info']['total_steps'],
                'best_mutant': workflow_data['results']['best_mutant']['name'],
                'improvement': workflow_data['results']['report']['assay_results']['improvement']
            })
        
        return {
            'status': 'success',
            'workflows': workflows,
            'total_workflows': len(workflows)
        }
    
    def handle_export_workflow_data(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle export workflow data command."""
        workflow_id = params.get('workflow_id')
        export_format = params.get('format', 'json')
        
        if not workflow_id:
            return {
                'status': 'error',
                'message': 'No workflow ID provided'
            }
        
        if workflow_id not in self.active_workflows:
            return {
                'status': 'error',
                'message': f'Workflow {workflow_id} not found'
            }
        
        workflow_data = self.active_workflows[workflow_id]
        
        if export_format == 'json':
            export_data = json.dumps(workflow_data, indent=2)
        else:
            export_data = str(workflow_data)
        
        return {
            'status': 'success',
            'workflow_id': workflow_id,
            'export_format': export_format,
            'export_data': export_data
        }


def register_directed_evolution_commands(handler: DirectedEvolutionHandler) -> Dict[str, Any]:
    """Register directed evolution commands with the command system."""
    commands = {
        'start_evolution': {
            'description': 'Start a directed evolution workflow',
            'parameters': {
                'target_property': 'Property to optimize (default: thermal_stability)',
                'library_size': 'Size of mutant library (default: 50)',
                'strategy': 'Mutation strategy (default: rational)',
                'assay_type': 'Type of assay to run (default: thermal_stability)',
                'name': 'Protein name (default: Directed_Evolution_Protein)',
                'description': 'Protein description',
                'length': 'Sequence length for generation (default: 250)'
            },
            'handler': handler.handle_start_evolution
        },
        'import_sequence_evolution': {
            'description': 'Run evolution with imported protein sequence',
            'parameters': {
                'sequence': 'Protein sequence to import',
                'target_property': 'Property to optimize (default: thermal_stability)',
                'library_size': 'Size of mutant library (default: 50)',
                'strategy': 'Mutation strategy (default: rational)',
                'assay_type': 'Type of assay to run (default: thermal_stability)',
                'name': 'Protein name (default: Imported_Protein)',
                'description': 'Protein description'
            },
            'handler': handler.handle_import_sequence_evolution
        },
        'multi_cycle_evolution': {
            'description': 'Run multi-cycle directed evolution',
            'parameters': {
                'num_cycles': 'Number of evolution cycles (default: 3)',
                'target_property': 'Property to optimize (default: thermal_stability)',
                'library_size': 'Initial library size (default: 30)',
                'strategy': 'Mutation strategy (default: rational)',
                'assay_type': 'Type of assay to run (default: thermal_stability)'
            },
            'handler': handler.handle_multi_cycle_evolution
        },
        'get_workflow_status': {
            'description': 'Get status of a directed evolution workflow',
            'parameters': {
                'workflow_id': 'ID of the workflow to check'
            },
            'handler': handler.handle_get_workflow_status
        },
        'list_workflows': {
            'description': 'List all active directed evolution workflows',
            'parameters': {},
            'handler': handler.handle_list_workflows
        },
        'export_workflow_data': {
            'description': 'Export workflow data',
            'parameters': {
                'workflow_id': 'ID of the workflow to export',
                'format': 'Export format (default: json)'
            },
            'handler': handler.handle_export_workflow_data
        }
    }
    
    return commands


if __name__ == "__main__":
    # Test the handler
    handler = DirectedEvolutionHandler()
    
    # Test start evolution
    print("Testing directed evolution handler...")
    
    result = handler.handle_start_evolution({
        'target_property': 'thermal_stability',
        'library_size': 30,
        'strategy': 'rational',
        'assay_type': 'thermal_stability'
    })
    
    print(f"Result: {result}") 