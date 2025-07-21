"""
Directed Evolution Workflow Engine Integration

This module provides workflow orchestration for directed evolution processes,
integrating with the platform's workflow engine and command system.
"""

import json
import uuid
from datetime import datetime
from typing import Dict, List, Any, Optional
import sys
import os

# Add the tools directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tools'))

from directed_evolution import DirectedEvolutionTool, ProteinSequence, EvolutionCycle


class DirectedEvolutionWorkflow:
    """Workflow orchestrator for directed evolution processes."""
    
    def __init__(self):
        self.tool = DirectedEvolutionTool()
        self.workflow_id = f"de_workflow_{uuid.uuid4().hex[:8]}"
        self.steps = []
        self.results = {}
        
    def start_workflow(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Start a directed evolution workflow."""
        workflow_info = {
            'workflow_id': self.workflow_id,
            'start_time': datetime.now().isoformat(),
            'config': config,
            'status': 'running'
        }
        
        self.results['workflow_info'] = workflow_info
        return workflow_info
    
    def step_1_generate_sequence(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 1: Generate or import protein sequence."""
        step_id = f"step_1_{uuid.uuid4().hex[:8]}"
        
        if config.get('import_sequence'):
            # Import existing sequence
            protein = self.tool.import_protein_sequence(
                sequence=config['sequence'],
                name=config.get('name', 'Imported_Protein'),
                description=config.get('description', 'Imported protein sequence')
            )
        else:
            # Generate new sequence
            protein = self.tool.generate_protein_sequence(
                name=config.get('name', 'Generated_Protein'),
                description=config.get('description', 'Generated protein sequence'),
                length=config.get('length', 300)
            )
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Generate/Import Sequence',
            'status': 'completed',
            'protein': {
                'id': protein.id,
                'name': protein.name,
                'length': len(protein.sequence),
                'properties': protein.properties
            },
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['sequence'] = protein
        
        return step_result
    
    def step_2_predict_mutations(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 2: Predict mutation hotspots."""
        step_id = f"step_2_{uuid.uuid4().hex[:8]}"
        
        protein = self.results['sequence']
        target_property = config.get('target_property', 'thermal_stability')
        
        mutations = self.tool.predict_mutation_hotspots(protein, target_property)
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Predict Mutations',
            'status': 'completed',
            'mutations': [
                {
                    'position': m.position,
                    'original_aa': m.original_aa,
                    'mutant_aa': m.mutant_aa,
                    'confidence': m.confidence,
                    'rationale': m.rationale
                } for m in mutations
            ],
            'target_property': target_property,
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['mutations'] = mutations
        
        return step_result
    
    def step_3_generate_library(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 3: Generate mutant library."""
        step_id = f"step_3_{uuid.uuid4().hex[:8]}"
        
        protein = self.results['sequence']
        strategy = config.get('strategy', 'rational')
        library_size = config.get('library_size', 50)
        
        library = self.tool.simulate_mutant_library(protein, strategy, library_size)
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Generate Library',
            'status': 'completed',
            'library': {
                'id': library.id,
                'name': library.name,
                'size': library.library_size,
                'strategy': library.strategy,
                'mutants': [
                    {
                        'id': m.id,
                        'name': m.name,
                        'properties': m.properties
                    } for m in library.mutants
                ]
            },
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['library'] = library
        
        return step_result
    
    def step_4_simulate_synthesis(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 4: Simulate cloning and synthesis."""
        step_id = f"step_4_{uuid.uuid4().hex[:8]}"
        
        library = self.results['library']
        synthesis_results = self.tool.simulate_cloning_synthesis(library)
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Simulate Synthesis',
            'status': 'completed',
            'synthesis_results': synthesis_results,
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['synthesis'] = synthesis_results
        
        return step_result
    
    def step_5_import_assay_data(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 5: Import assay data."""
        step_id = f"step_5_{uuid.uuid4().hex[:8]}"
        
        library = self.results['library']
        assay_type = config.get('assay_type', 'thermal_stability')
        
        assay_results = self.tool.import_assay_data(library, assay_type)
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Import Assay Data',
            'status': 'completed',
            'assay_results': [
                {
                    'mutant_id': r.mutant_id,
                    'assay_type': r.assay_type,
                    'value': r.value,
                    'unit': r.unit,
                    'confidence': r.confidence
                } for r in assay_results
            ],
            'assay_type': assay_type,
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['assay_results'] = assay_results
        
        return step_result
    
    def step_6_analyze_results(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 6: Analyze sequence-phenotype relationships."""
        step_id = f"step_6_{uuid.uuid4().hex[:8]}"
        
        library = self.results['library']
        assay_results = self.results['assay_results']
        
        analysis = self.tool.analyze_sequence_phenotype_relationships(library, assay_results)
        best_mutant = self.tool.identify_best_mutant(library, assay_results)
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Analyze Results',
            'status': 'completed',
            'analysis': analysis,
            'best_mutant': {
                'id': best_mutant.id,
                'name': best_mutant.name,
                'properties': best_mutant.properties
            },
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['analysis'] = analysis
        self.results['best_mutant'] = best_mutant
        
        return step_result
    
    def step_7_generate_report(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Step 7: Generate comprehensive report."""
        step_id = f"step_7_{uuid.uuid4().hex[:8]}"
        
        # Create a cycle object for reporting
        cycle = EvolutionCycle(
            cycle_number=1,
            input_sequence=self.results['sequence'],
            generated_library=self.results['library'],
            assay_results=self.results['assay_results'],
            best_mutant=self.results['best_mutant'],
            improvement_metrics={
                'property': config.get('target_property', 'thermal_stability'),
                'original_value': self.results['sequence'].properties.get('thermal_stability', 0),
                'best_value': max(r.value for r in self.results['assay_results']),
                'improvement': 0,  # Will be calculated
                'improvement_percentage': 0  # Will be calculated
            },
            learnings=[],
            next_steps=[],
            created_at=datetime.now().isoformat()
        )
        
        # Calculate improvement
        original_value = cycle.improvement_metrics['original_value']
        best_value = cycle.improvement_metrics['best_value']
        improvement = best_value - original_value
        cycle.improvement_metrics['improvement'] = improvement
        cycle.improvement_metrics['improvement_percentage'] = (improvement / original_value * 100) if original_value > 0 else 0
        
        report = self.tool.generate_evolution_report(cycle)
        visualization = self.tool.create_visualization(cycle)
        export_data = self.tool.export_cycle_data(cycle)
        
        step_result = {
            'step_id': step_id,
            'step_name': 'Generate Report',
            'status': 'completed',
            'report': report,
            'visualization': visualization,
            'export_data': export_data,
            'timestamp': datetime.now().isoformat()
        }
        
        self.steps.append(step_result)
        self.results['report'] = report
        self.results['visualization'] = visualization
        self.results['export_data'] = export_data
        
        return step_result
    
    def run_complete_workflow(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Run the complete directed evolution workflow."""
        print("🧬 Starting Directed Evolution Workflow")
        print("=" * 50)
        
        # Start workflow
        workflow_info = self.start_workflow(config)
        print(f"✅ Workflow started: {workflow_info['workflow_id']}")
        
        # Step 1: Generate/Import sequence
        print("\n1. Generating/Importing sequence...")
        step1 = self.step_1_generate_sequence(config)
        print(f"✅ {step1['step_name']} completed")
        print(f"   Protein: {step1['protein']['name']} ({step1['protein']['length']} aa)")
        
        # Step 2: Predict mutations
        print("\n2. Predicting mutations...")
        step2 = self.step_2_predict_mutations(config)
        print(f"✅ {step2['step_name']} completed")
        print(f"   Mutations predicted: {len(step2['mutations'])}")
        
        # Step 3: Generate library
        print("\n3. Generating mutant library...")
        step3 = self.step_3_generate_library(config)
        print(f"✅ {step3['step_name']} completed")
        print(f"   Library size: {step3['library']['size']}")
        
        # Step 4: Simulate synthesis
        print("\n4. Simulating synthesis...")
        step4 = self.step_4_simulate_synthesis(config)
        print(f"✅ {step4['step_name']} completed")
        print(f"   Successful clones: {step4['synthesis_results']['successful_clones']}")
        
        # Step 5: Import assay data
        print("\n5. Importing assay data...")
        step5 = self.step_5_import_assay_data(config)
        print(f"✅ {step5['step_name']} completed")
        print(f"   Assay results: {len(step5['assay_results'])}")
        
        # Step 6: Analyze results
        print("\n6. Analyzing results...")
        step6 = self.step_6_analyze_results(config)
        print(f"✅ {step6['step_name']} completed")
        print(f"   Best mutant: {step6['best_mutant']['name']}")
        
        # Step 7: Generate report
        print("\n7. Generating report...")
        step7 = self.step_7_generate_report(config)
        print(f"✅ {step7['step_name']} completed")
        
        # Complete workflow
        workflow_info['status'] = 'completed'
        workflow_info['end_time'] = datetime.now().isoformat()
        workflow_info['total_steps'] = len(self.steps)
        
        print("\n" + "=" * 50)
        print("🎉 Directed Evolution Workflow Complete!")
        print(f"   Workflow ID: {workflow_info['workflow_id']}")
        print(f"   Total steps: {workflow_info['total_steps']}")
        print(f"   Best improvement: {step7['report']['assay_results']['improvement']:.2f}")
        
        return {
            'workflow_info': workflow_info,
            'steps': self.steps,
            'results': self.results
        }


def create_workflow_config(target_property: str = "thermal_stability",
                          library_size: int = 50,
                          strategy: str = "rational",
                          assay_type: str = "thermal_stability") -> Dict[str, Any]:
    """Create a standard workflow configuration."""
    return {
        'target_property': target_property,
        'library_size': library_size,
        'strategy': strategy,
        'assay_type': assay_type,
        'name': 'Directed_Evolution_Protein',
        'description': 'Protein for directed evolution optimization',
        'length': 250,
        'import_sequence': False
    }


if __name__ == "__main__":
    # Example usage
    config = create_workflow_config()
    workflow = DirectedEvolutionWorkflow()
    results = workflow.run_complete_workflow(config)
    
    print("\nWorkflow Results Summary:")
    print(f"Best mutant: {results['results']['best_mutant']['name']}")
    print(f"Improvement: {results['results']['report']['assay_results']['improvement']:.2f}")
    print(f"Total mutants: {results['results']['library']['size']}") 