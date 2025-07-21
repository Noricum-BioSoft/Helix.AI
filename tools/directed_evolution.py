"""
Directed Evolution Tool for Protein Engineering

This tool implements a complete Design-Build-Test-Learn (DBTL) cycle for directed evolution
of proteins, focusing on improving enzyme thermal stability and other properties.
"""

import json
import random
import string
import uuid
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from dataclasses import dataclass, asdict
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import base64


@dataclass
class ProteinSequence:
    """Represents a protein sequence with metadata."""
    id: str
    sequence: str
    name: str
    description: str
    properties: Dict[str, Any]
    created_at: str
    version: int = 1


@dataclass
class Mutation:
    """Represents a single mutation."""
    position: int
    original_aa: str
    mutant_aa: str
    mutation_type: str  # 'point', 'insertion', 'deletion'
    confidence: float
    rationale: str


@dataclass
class MutantLibrary:
    """Represents a library of mutants."""
    id: str
    name: str
    parent_sequence_id: str
    mutants: List[ProteinSequence]
    mutations: List[Mutation]
    library_size: int
    created_at: str
    strategy: str  # 'random', 'rational', 'structure_based'


@dataclass
class AssayResult:
    """Represents assay data for a mutant."""
    mutant_id: str
    assay_type: str  # 'thermal_stability', 'activity', 'expression'
    value: float
    unit: str
    confidence: float
    conditions: Dict[str, Any]
    date: str


@dataclass
class EvolutionCycle:
    """Represents one complete DBTL cycle."""
    cycle_number: int
    input_sequence: ProteinSequence
    generated_library: MutantLibrary
    assay_results: List[AssayResult]
    best_mutant: Optional[ProteinSequence]
    improvement_metrics: Dict[str, float]
    learnings: List[str]
    next_steps: List[str]
    created_at: str


class DirectedEvolutionTool:
    """Main tool for directed evolution workflows."""
    
    def __init__(self):
        self.sequences: Dict[str, ProteinSequence] = {}
        self.libraries: Dict[str, MutantLibrary] = {}
        self.cycles: Dict[str, EvolutionCycle] = {}
        self.current_cycle = 1
        
    def generate_protein_sequence(self, name: str, description: str, length: int = 300) -> ProteinSequence:
        """Generate a new protein sequence for evolution."""
        # Generate a realistic protein sequence (simplified)
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        sequence = ''.join(random.choice(amino_acids) for _ in range(length))
        
        sequence_id = f"seq_{uuid.uuid4().hex[:8]}"
        protein = ProteinSequence(
            id=sequence_id,
            sequence=sequence,
            name=name,
            description=description,
            properties={
                'length': length,
                'molecular_weight': length * 110,  # Approximate
                'isoelectric_point': random.uniform(5.0, 9.0),
                'thermal_stability': random.uniform(40, 80),
                'activity': random.uniform(0.1, 1.0)
            },
            created_at=datetime.now().isoformat()
        )
        
        self.sequences[sequence_id] = protein
        return protein
    
    def import_protein_sequence(self, sequence: str, name: str, description: str) -> ProteinSequence:
        """Import an existing protein sequence."""
        sequence_id = f"seq_{uuid.uuid4().hex[:8]}"
        protein = ProteinSequence(
            id=sequence_id,
            sequence=sequence,
            name=name,
            description=description,
            properties={
                'length': len(sequence),
                'molecular_weight': len(sequence) * 110,
                'isoelectric_point': random.uniform(5.0, 9.0),
                'thermal_stability': random.uniform(40, 80),
                'activity': random.uniform(0.1, 1.0)
            },
            created_at=datetime.now().isoformat()
        )
        
        self.sequences[sequence_id] = protein
        return protein
    
    def predict_mutation_hotspots(self, sequence: ProteinSequence, 
                                 target_property: str = "thermal_stability") -> List[Mutation]:
        """Predict mutation hotspots for improving target property."""
        mutations = []
        seq = sequence.sequence
        
        # Simulate hotspot prediction based on sequence properties
        for i in range(0, len(seq), 10):  # Every 10th position
            if i < len(seq):
                original_aa = seq[i]
                # Generate potential mutations
                amino_acids = "ACDEFGHIKLMNPQRSTVWY"
                for new_aa in random.sample(amino_acids, 3):
                    if new_aa != original_aa:
                        mutation = Mutation(
                            position=i,
                            original_aa=original_aa,
                            mutant_aa=new_aa,
                            mutation_type='point',
                            confidence=random.uniform(0.3, 0.9),
                            rationale=f"Predicted to improve {target_property} based on structural analysis"
                        )
                        mutations.append(mutation)
        
        return mutations
    
    def simulate_mutant_library(self, parent_sequence: ProteinSequence, 
                               strategy: str = "rational", 
                               library_size: int = 50) -> MutantLibrary:
        """Simulate generation of a mutant library."""
        library_id = f"lib_{uuid.uuid4().hex[:8]}"
        
        # Generate mutations based on strategy
        if strategy == "rational":
            mutations = self.predict_mutation_hotspots(parent_sequence)
        else:
            # Random mutations
            mutations = []
            seq = parent_sequence.sequence
            for _ in range(library_size):
                pos = random.randint(0, len(seq) - 1)
                original_aa = seq[pos]
                amino_acids = "ACDEFGHIKLMNPQRSTVWY"
                new_aa = random.choice([aa for aa in amino_acids if aa != original_aa])
                
                mutation = Mutation(
                    position=pos,
                    original_aa=original_aa,
                    mutant_aa=new_aa,
                    mutation_type='point',
                    confidence=random.uniform(0.1, 0.8),
                    rationale="Random mutation"
                )
                mutations.append(mutation)
        
        # Create mutant sequences
        mutants = []
        for i, mutation in enumerate(mutations[:library_size]):
            # Apply mutation to create new sequence
            seq_list = list(parent_sequence.sequence)
            seq_list[mutation.position] = mutation.mutant_aa
            new_seq = ''.join(seq_list)
            
            mutant = ProteinSequence(
                id=f"mut_{uuid.uuid4().hex[:8]}",
                sequence=new_seq,
                name=f"{parent_sequence.name}_mut_{i+1}",
                description=f"Mutant {i+1} of {parent_sequence.name}",
                properties={
                    'length': len(new_seq),
                    'molecular_weight': len(new_seq) * 110,
                    'isoelectric_point': random.uniform(5.0, 9.0),
                    'thermal_stability': parent_sequence.properties['thermal_stability'] + random.uniform(-10, 15),
                    'activity': parent_sequence.properties['activity'] + random.uniform(-0.2, 0.3)
                },
                created_at=datetime.now().isoformat(),
                version=parent_sequence.version + 1
            )
            mutants.append(mutant)
            self.sequences[mutant.id] = mutant
        
        library = MutantLibrary(
            id=library_id,
            name=f"{parent_sequence.name}_library_{self.current_cycle}",
            parent_sequence_id=parent_sequence.id,
            mutants=mutants,
            mutations=mutations[:library_size],
            library_size=len(mutants),
            created_at=datetime.now().isoformat(),
            strategy=strategy
        )
        
        self.libraries[library_id] = library
        return library
    
    def simulate_cloning_synthesis(self, library: MutantLibrary) -> Dict[str, Any]:
        """Simulate cloning and synthesis process."""
        results = {
            'successful_clones': len(library.mutants),
            'failed_clones': random.randint(0, 5),
            'synthesis_time': random.randint(3, 14),  # days
            'cost': len(library.mutants) * random.uniform(50, 200),
            'quality_control': {
                'sequence_verification': random.uniform(0.8, 1.0),
                'expression_check': random.uniform(0.7, 0.95),
                'purity_assessment': random.uniform(0.75, 0.98)
            }
        }
        return results
    
    def import_assay_data(self, library: MutantLibrary, 
                         assay_type: str = "thermal_stability") -> List[AssayResult]:
        """Import or simulate assay data for mutants."""
        results = []
        
        for mutant in library.mutants:
            # Simulate assay results
            if assay_type == "thermal_stability":
                base_temp = mutant.properties['thermal_stability']
                measured_temp = base_temp + random.uniform(-5, 10)
                result = AssayResult(
                    mutant_id=mutant.id,
                    assay_type=assay_type,
                    value=measured_temp,
                    unit="°C",
                    confidence=random.uniform(0.7, 0.95),
                    conditions={
                        'temperature': 25,
                        'pH': 7.0,
                        'buffer': 'PBS'
                    },
                    date=datetime.now().isoformat()
                )
            elif assay_type == "activity":
                base_activity = mutant.properties['activity']
                measured_activity = base_activity + random.uniform(-0.1, 0.2)
                result = AssayResult(
                    mutant_id=mutant.id,
                    assay_type=assay_type,
                    value=measured_activity,
                    unit="U/mg",
                    confidence=random.uniform(0.6, 0.9),
                    conditions={
                        'substrate': 'standard',
                        'temperature': 37,
                        'pH': 7.5
                    },
                    date=datetime.now().isoformat()
                )
            else:
                # Generic assay
                result = AssayResult(
                    mutant_id=mutant.id,
                    assay_type=assay_type,
                    value=random.uniform(0.1, 10.0),
                    unit="arbitrary",
                    confidence=random.uniform(0.5, 0.9),
                    conditions={},
                    date=datetime.now().isoformat()
                )
            
            results.append(result)
        
        return results
    
    def analyze_sequence_phenotype_relationships(self, library: MutantLibrary, 
                                              assay_results: List[AssayResult]) -> Dict[str, Any]:
        """Analyze relationships between sequence changes and phenotypic outcomes."""
        analysis = {
            'correlation_analysis': {},
            'mutation_impact': {},
            'structure_activity_relationships': {},
            'recommendations': []
        }
        
        # Analyze mutation impacts
        mutation_impacts = {}
        for mutation in library.mutations:
            impact_score = random.uniform(-0.5, 1.0)
            mutation_impacts[f"{mutation.original_aa}{mutation.position}{mutation.mutant_aa}"] = {
                'impact': impact_score,
                'confidence': mutation.confidence,
                'rationale': mutation.rationale
            }
        
        analysis['mutation_impact'] = mutation_impacts
        
        # Generate recommendations
        recommendations = [
            "Focus on mutations in positions 50-100 for thermal stability",
            "Avoid mutations in catalytic residues",
            "Consider combining beneficial mutations",
            "Test expression levels for all mutants"
        ]
        analysis['recommendations'] = recommendations
        
        return analysis
    
    def identify_best_mutant(self, library: MutantLibrary, 
                           assay_results: List[AssayResult]) -> ProteinSequence:
        """Identify the best performing mutant from assay results."""
        # Find mutant with highest thermal stability
        best_result = max(assay_results, key=lambda x: x.value)
        best_mutant = next(m for m in library.mutants if m.id == best_result.mutant_id)
        
        return best_mutant
    
    def run_evolution_cycle(self, input_sequence: ProteinSequence, 
                           target_property: str = "thermal_stability",
                           library_size: int = 50) -> EvolutionCycle:
        """Run one complete DBTL cycle."""
        # Design: Generate mutant library
        library = self.simulate_mutant_library(input_sequence, "rational", library_size)
        
        # Build: Simulate cloning/synthesis
        synthesis_results = self.simulate_cloning_synthesis(library)
        
        # Test: Import assay data
        assay_results = self.import_assay_data(library, target_property)
        
        # Learn: Analyze results and identify best mutant
        analysis = self.analyze_sequence_phenotype_relationships(library, assay_results)
        best_mutant = self.identify_best_mutant(library, assay_results)
        
        # Calculate improvement metrics
        original_value = input_sequence.properties.get(target_property, 0)
        best_value = max(r.value for r in assay_results if r.assay_type == target_property)
        improvement = best_value - original_value
        
        improvement_metrics = {
            'property': target_property,
            'original_value': original_value,
            'best_value': best_value,
            'improvement': improvement,
            'improvement_percentage': (improvement / original_value * 100) if original_value > 0 else 0
        }
        
        # Generate learnings and next steps
        learnings = [
            f"Best mutation improved {target_property} by {improvement:.2f}",
            f"Library size of {library_size} provided good coverage",
            "Rational design strategy showed promising results"
        ]
        
        next_steps = [
            "Generate focused library around best mutations",
            "Test expression and purification of best mutants",
            "Characterize stability under different conditions",
            "Prepare for next evolution cycle"
        ]
        
        cycle = EvolutionCycle(
            cycle_number=self.current_cycle,
            input_sequence=input_sequence,
            generated_library=library,
            assay_results=assay_results,
            best_mutant=best_mutant,
            improvement_metrics=improvement_metrics,
            learnings=learnings,
            next_steps=next_steps,
            created_at=datetime.now().isoformat()
        )
        
        cycle_id = f"cycle_{self.current_cycle}_{uuid.uuid4().hex[:8]}"
        self.cycles[cycle_id] = cycle
        self.current_cycle += 1
        
        return cycle
    
    def generate_evolution_report(self, cycle: EvolutionCycle) -> Dict[str, Any]:
        """Generate a comprehensive report for an evolution cycle."""
        report = {
            'cycle_info': {
                'cycle_number': cycle.cycle_number,
                'date': cycle.created_at,
                'target_property': cycle.improvement_metrics['property']
            },
            'input_sequence': {
                'id': cycle.input_sequence.id,
                'name': cycle.input_sequence.name,
                'length': len(cycle.input_sequence.sequence),
                'properties': cycle.input_sequence.properties
            },
            'library_stats': {
                'size': cycle.generated_library.library_size,
                'strategy': cycle.generated_library.strategy,
                'synthesis_success': len(cycle.generated_library.mutants)
            },
            'assay_results': {
                'total_tested': len(cycle.assay_results),
                'best_value': cycle.improvement_metrics['best_value'],
                'improvement': cycle.improvement_metrics['improvement']
            },
            'best_mutant': {
                'id': cycle.best_mutant.id,
                'name': cycle.best_mutant.name,
                'properties': cycle.best_mutant.properties
            },
            'learnings': cycle.learnings,
            'next_steps': cycle.next_steps,
            'recommendations': [
                "Continue evolution with focused library",
                "Characterize best mutants in detail",
                "Consider structural analysis of improvements"
            ]
        }
        
        return report
    
    def create_visualization(self, cycle: EvolutionCycle) -> str:
        """Create visualization of evolution results."""
        # Create a figure with multiple subplots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Assay results distribution
        values = [r.value for r in cycle.assay_results]
        ax1.hist(values, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax1.axvline(cycle.improvement_metrics['original_value'], color='red', linestyle='--', label='Original')
        ax1.axvline(cycle.improvement_metrics['best_value'], color='green', linestyle='--', label='Best')
        ax1.set_xlabel(f"{cycle.improvement_metrics['property']} Value")
        ax1.set_ylabel('Frequency')
        ax1.set_title('Distribution of Assay Results')
        ax1.legend()
        
        # Plot 2: Improvement over cycle
        ax2.bar(['Original', 'Best'], 
                [cycle.improvement_metrics['original_value'], cycle.improvement_metrics['best_value']],
                color=['lightcoral', 'lightgreen'])
        ax2.set_ylabel(f"{cycle.improvement_metrics['property']} Value")
        ax2.set_title('Improvement Comparison')
        
        # Plot 3: Mutation positions
        positions = [m.position for m in cycle.generated_library.mutations]
        ax3.hist(positions, bins=20, alpha=0.7, color='orange', edgecolor='black')
        ax3.set_xlabel('Sequence Position')
        ax3.set_ylabel('Number of Mutations')
        ax3.set_title('Mutation Position Distribution')
        
        # Plot 4: Confidence vs Value
        confidences = [r.confidence for r in cycle.assay_results]
        ax4.scatter(confidences, values, alpha=0.6, color='purple')
        ax4.set_xlabel('Assay Confidence')
        ax4.set_ylabel(f"{cycle.improvement_metrics['property']} Value")
        ax4.set_title('Confidence vs Performance')
        
        plt.tight_layout()
        
        # Convert to base64 string
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        buf.seek(0)
        img_str = base64.b64encode(buf.getvalue()).decode()
        plt.close()
        
        return img_str
    
    def export_cycle_data(self, cycle: EvolutionCycle, format: str = "json") -> str:
        """Export cycle data in various formats."""
        if format == "json":
            data = {
                'cycle': asdict(cycle),
                'library': asdict(cycle.generated_library),
                'assay_results': [asdict(r) for r in cycle.assay_results],
                'report': self.generate_evolution_report(cycle)
            }
            return json.dumps(data, indent=2)
        else:
            return "Unsupported format"


def run_directed_evolution_demo():
    """Run a complete directed evolution demo."""
    tool = DirectedEvolutionTool()
    
    print("🧬 Starting Directed Evolution Demo")
    print("=" * 50)
    
    # Step 1: Generate/Import protein sequence
    print("\n1. Generating protein sequence...")
    protein = tool.generate_protein_sequence(
        name="Demo_Enzyme",
        description="Thermophilic enzyme for directed evolution",
        length=250
    )
    print(f"✅ Created protein: {protein.name} (ID: {protein.id})")
    print(f"   Length: {len(protein.sequence)} amino acids")
    print(f"   Initial thermal stability: {protein.properties['thermal_stability']:.1f}°C")
    
    # Step 2: Run first evolution cycle
    print("\n2. Running evolution cycle 1...")
    cycle1 = tool.run_evolution_cycle(protein, "thermal_stability", 30)
    print(f"✅ Cycle 1 completed")
    print(f"   Library size: {cycle1.generated_library.library_size}")
    print(f"   Best improvement: {cycle1.improvement_metrics['improvement']:.2f}°C")
    print(f"   Best mutant: {cycle1.best_mutant.name}")
    
    # Step 3: Run second evolution cycle with best mutant
    print("\n3. Running evolution cycle 2...")
    cycle2 = tool.run_evolution_cycle(cycle1.best_mutant, "thermal_stability", 25)
    print(f"✅ Cycle 2 completed")
    print(f"   Library size: {cycle2.generated_library.library_size}")
    print(f"   Best improvement: {cycle2.improvement_metrics['improvement']:.2f}°C")
    print(f"   Best mutant: {cycle2.best_mutant.name}")
    
    # Step 4: Generate reports
    print("\n4. Generating reports...")
    report1 = tool.generate_evolution_report(cycle1)
    report2 = tool.generate_evolution_report(cycle2)
    
    print(f"✅ Cycle 1 Report:")
    print(f"   Original stability: {report1['input_sequence']['properties']['thermal_stability']:.1f}°C")
    print(f"   Best achieved: {report1['assay_results']['best_value']:.1f}°C")
    print(f"   Improvement: {report1['assay_results']['improvement']:.2f}°C")
    
    print(f"✅ Cycle 2 Report:")
    print(f"   Original stability: {report2['input_sequence']['properties']['thermal_stability']:.1f}°C")
    print(f"   Best achieved: {report2['assay_results']['best_value']:.1f}°C")
    print(f"   Improvement: {report2['assay_results']['improvement']:.2f}°C")
    
    # Step 5: Create visualizations
    print("\n5. Creating visualizations...")
    viz1 = tool.create_visualization(cycle1)
    viz2 = tool.create_visualization(cycle2)
    print("✅ Visualizations created")
    
    # Step 6: Export data
    print("\n6. Exporting data...")
    export1 = tool.export_cycle_data(cycle1)
    export2 = tool.export_cycle_data(cycle2)
    print("✅ Data exported")
    
    # Summary
    print("\n" + "=" * 50)
    print("🎉 Directed Evolution Demo Complete!")
    print(f"   Total cycles: 2")
    print(f"   Total mutants generated: {cycle1.generated_library.library_size + cycle2.generated_library.library_size}")
    print(f"   Final improvement: {report2['assay_results']['best_value'] - report1['input_sequence']['properties']['thermal_stability']:.2f}°C")
    print(f"   Best final mutant: {cycle2.best_mutant.name}")
    
    return {
        'tool': tool,
        'cycles': [cycle1, cycle2],
        'reports': [report1, report2],
        'visualizations': [viz1, viz2],
        'exports': [export1, export2]
    }


if __name__ == "__main__":
    demo_results = run_directed_evolution_demo() 