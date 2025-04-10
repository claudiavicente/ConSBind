#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Add these imports at the top of your file
from Bio.PDB import PDBParser, Selection
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import logging

def pred_scores(prediction_file):
    predicted_residues = []
    scores = {
        'consensus_score': None,
        'binding_potential_score': None,
        'druggability_score': None,
        'knowledge_based_score': None
    }
    
    with open(prediction_file, 'r') as f:
        lines = f.readlines()
        
        in_residues_section = False
        for line in lines:
            line = line.strip()
            
            # Extract scores
            if line.startswith("Consensus Score:"):
                scores['consensus_score'] = float(line.split(":")[1].strip())
            elif line.startswith("Binding Potential Score:"):
                scores['binding_potential_score'] = float(line.split(":")[1].strip())
            elif line.startswith("Druggability:"):
                scores['druggability_score'] = float(line.split(":")[1].strip())
            elif line.startswith("Knowledge-based Score:"):
                scores['knowledge_based_score'] = float(line.split(":")[1].strip())
            
            if line == "Binding Site Residues:":
                in_residues_section = True
                continue
            
            if in_residues_section and line and not line.startswith("Site"):
                if line == "":
                    in_residues_section = False
                    continue
                
                parts = line.strip().split(':')
                if len(parts) == 2:
                    chain_id = parts[0].strip()
                    residue_info = parts[1].strip()
                    
                    residue_id = ''.join(filter(str.isdigit, residue_info))
                    if residue_id:
                        predicted_residues.append((chain_id, int(residue_id)))
    
    return predicted_residues, scores

def calculate_complexity(pdb_id, structure, structure_file, residue_map, known_residues):
    complexity_score = 0.0
    
    try:
        # 1. Size complexity - larger proteins are more complex
        all_residues = Selection.unfold_entities(structure, 'R')
        standard_aa_names = set(["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", 
                                "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", 
                                "SER", "THR", "VAL", "TRP", "TYR"])
        residue_count = len([r for r in all_residues if r.get_resname() in standard_aa_names])
        size_factor = min(1.0, residue_count / 500)  # Normalize by 500 residues
        complexity_score += 0.25 * size_factor
        
        # 2. Domain complexity - multi-domain proteins are more complex
        chains = list(structure.get_chains())
        chain_count = len(chains)
        chain_factor = min(1.0, chain_count / 4)  # Normalize by 4 chains
        complexity_score += 0.25 * chain_factor
        
        # 3. Structural complexity based on B-factors (indicates flexibility)
        b_factors = []
        for atom in structure.get_atoms():
            if atom.get_name() == 'CA':  # Only alpha carbons
                b_factors.append(atom.get_bfactor())
        
        if b_factors:
            # Higher B-factor variance indicates more complex/flexible structure
            b_factor_std = np.std(b_factors)
            b_factor_complexity = min(1.0, b_factor_std / 30)  # Normalize
            complexity_score += 0.25 * b_factor_complexity
        
        # 4. Binding site complexity
        # Calculate how buried or surface-exposed the binding site is
        if known_residues:
            # Get all atoms in the structure for distance calculations
            all_atoms = list(structure.get_atoms())
            all_coords = np.array([atom.get_coord() for atom in all_atoms])
            
            # Get binding site residue atoms
            binding_site_atoms = []
            for chain, res_id in known_residues:
                if (chain, res_id) in residue_map:
                    binding_site_atoms.extend(list(residue_map[(chain, res_id)].get_atoms()))
            
            if binding_site_atoms:
                # Calculate how buried each binding site atom is
                buried_scores = []
                for bs_atom in binding_site_atoms:
                    bs_coord = bs_atom.get_coord()
                    
                    # Count how many atoms are within 10Å
                    distances = np.linalg.norm(all_coords - bs_coord, axis=1)
                    nearby_atoms = sum(distances < 10.0)
                    
                    # Normalize by a typical fully buried value
                    buried_score = min(1.0, nearby_atoms / 300)
                    buried_scores.append(buried_score)
                
                # Average burial score for binding site
                avg_burial = np.mean(buried_scores) if buried_scores else 0.0
                complexity_score += 0.25 * avg_burial
        
        # Normalize to 0-10 scale
        complexity_score = min(1.0, complexity_score) * 10
        
    except Exception as e:
        logger.error(f"Error calculating complexity for {pdb_id}: {e}")
        complexity_score = None
    
    return complexity_score
    
def calculate_metrics(protein_class, pdb_id, known_residues, predicted_residues, prediction_scores):
    """
    Calculate performance metrics for a single protein
    """
    # Convert residues to sets for easier comparison
    known_set = set((chain, res_id) for chain, res_id in known_residues)
    predicted_set = set((chain, res_id) for chain, res_id in predicted_residues)
    
    # Calculate precision, recall, and F1 score
    true_positives = len(known_set.intersection(predicted_set))
    false_positives = len(predicted_set - known_set)
    false_negatives = len(known_set - predicted_set)
    
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    # Load structure to calculate spatial metrics
    structure_file = base_dir / protein_class / f"pdb{pdb_id}.ent"
    
    # Initialize values for metrics that require structure analysis
    spatial_overlap = None
    center_distance = None
    volume_similarity = None
    pocket_depth = None
    pocket_polarity = None
    complexity_score = None
    
    if structure_file.exists():
        # Calculate spatial metrics using the structure
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, structure_file)
            
            # Get atoms for known and predicted binding sites
            known_atoms = []
            predicted_atoms = []
            all_residues = Selection.unfold_entities(structure, 'R')
            
            # Map residue IDs to actual residue objects
            residue_map = {}
            for residue in all_residues:
                if residue.get_resname() in standard_aa_names:  # Only consider standard amino acids
                    chain_id = residue.get_parent().id
                    res_id = residue.get_id()[1]
                    residue_map[(chain_id, res_id)] = residue
            
            # Get atoms for known binding site residues
            for res_id in known_residues:
                if res_id in residue_map:
                    for atom in residue_map[res_id]:
                        if atom.get_name() not in ['H', 'HA']:  # Skip hydrogen atoms
                            known_atoms.append(atom)
            
            # Get atoms for predicted binding site residues
            for res_id in predicted_residues:
                if res_id in residue_map:
                    for atom in residue_map[res_id]:
                        if atom.get_name() not in ['H', 'HA']:  # Skip hydrogen atoms
                            predicted_atoms.append(atom)
            
            if len(known_atoms) > 0 and len(predicted_atoms) > 0:
                # 1. Calculate spatial overlap (Jaccard index in 3D space)
                # Create a grid representation of the binding sites
                grid_spacing = 1.0  # Angstroms
                known_grid = set()
                predicted_grid = set()
                
                # Define a grid around each atom
                for atom in known_atoms:
                    coord = atom.get_coord()
                    grid_point = tuple(np.round(coord / grid_spacing).astype(int))
                    known_grid.add(grid_point)
                    
                    # Add neighboring grid points within van der Waals radius
                    radius = 1.5  # Approximate vdW radius in grid units
                    for dx in range(-int(radius), int(radius)+1):
                        for dy in range(-int(radius), int(radius)+1):
                            for dz in range(-int(radius), int(radius)+1):
                                if dx*dx + dy*dy + dz*dz <= radius*radius:
                                    known_grid.add((grid_point[0]+dx, grid_point[1]+dy, grid_point[2]+dz))
                
                for atom in predicted_atoms:
                    coord = atom.get_coord()
                    grid_point = tuple(np.round(coord / grid_spacing).astype(int))
                    predicted_grid.add(grid_point)
                    
                    # Add neighboring grid points within van der Waals radius
                    radius = 1.5  # Approximate vdW radius in grid units
                    for dx in range(-int(radius), int(radius)+1):
                        for dy in range(-int(radius), int(radius)+1):
                            for dz in range(-int(radius), int(radius)+1):
                                if dx*dx + dy*dy + dz*dz <= radius*radius:
                                    predicted_grid.add((grid_point[0]+dx, grid_point[1]+dy, grid_point[2]+dz))

                
                complexity_score = calculate_complexity(pdb_id, structure, structure_file, residue_map, known_residues)

                # Calculate Jaccard index (intersection over union)
                intersection = len(known_grid.intersection(predicted_grid))
                union = len(known_grid.union(predicted_grid))
                spatial_overlap = intersection / union if union > 0 else 0.0
                
                # Calculate center distance
                known_coords = np.array([atom.get_coord() for atom in known_atoms])
                predicted_coords = np.array([atom.get_coord() for atom in predicted_atoms])
                known_center = np.mean(known_coords, axis=0)
                predicted_center = np.mean(predicted_coords, axis=0)
                center_distance = np.linalg.norm(known_center - predicted_center)
                
                # Calculate volume similarity
                known_volume = len(known_grid) * (grid_spacing ** 3)
                predicted_volume = len(predicted_grid) * (grid_spacing ** 3)
                volume_similarity = min(known_volume, predicted_volume) / max(known_volume, predicted_volume)
                
                # Calculate pocket depth
                structure_file = str(convert_ent_to_pdb(structure_file))
                protein = ProteinStructure(str(structure_file))
                surface_atoms = protein.get_surface_atoms()
                surface_coords = np.array([atom.get_coord() for atom in surface_atoms])
                if len(surface_coords) > 0:
                    distances = np.linalg.norm(surface_coords[:, np.newaxis] - predicted_center, axis=2)
                    pocket_depth = np.min(distances)
                
                # Calculate pocket polarity
                polar_atoms = ['N', 'O', 'S']
                polar_count = sum(1 for atom in predicted_atoms if atom.get_name()[0] in polar_atoms)
                pocket_polarity = polar_count / len(predicted_atoms) if len(predicted_atoms) > 0 else 0.0
        except Exception as e:
            logger.error(f"Error calculating spatial metrics for {pdb_id}: {e}")
    
    return {
        'protein_complexity': complexity_score,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'spatial_overlap': spatial_overlap,
        'center_distance': center_distance,
        'volume_similarity': volume_similarity,
        'pocket_depth': pocket_depth,
        'pocket_polarity': pocket_polarity,
        'num_predictions': len(predicted_residues),
        'consensus_score': prediction_scores['consensus_score'],
        'binding_potential_score': prediction_scores['binding_potential_score'],
        'druggability_score': prediction_scores['druggability_score'],
        'knowledge_based_score': prediction_scores['knowledge_based_score']
    }

def evaluate_predictions():

    results = []
    
    for protein_class, proteins in PROTEIN_CLASSES.items():
        logger.info(f"Evaluating {len(proteins)} proteins for class: {protein_class}")
        
        for pdb_id in tqdm(proteins.keys(), desc=f"Evaluating {protein_class}"):
            # Find prediction results
            results_dir = Path('results/analysis')
            if not results_dir.exists():
                logger.warning(f"Results directory not found: {results_dir}")
                continue
            
            # Look for prediction files
            prediction_files = list(results_dir.glob(f"**/pdb{pdb_id}_predictions.txt"))
            if not prediction_files:
                logger.warning(f"No prediction file found for {pdb_id}")
                continue
            
            prediction_file = prediction_files[0]
            
            # Parse prediction file to get predicted binding sites and scores
            predicted_residues, prediction_scores = pred_scores(prediction_file)
            
            # Get known binding site residues
            known_residues = KNOWN_BINDING_SITES.get(pdb_id, [])
            
            if not known_residues:
                logger.warning(f"No known binding site residues for {pdb_id}")
                continue
            
            # Calculate metrics
            metrics = calculate_metrics(protein_class, pdb_id, known_residues, predicted_residues, prediction_scores)
            
            # Add to results
            results.append({
                'pdb_id': pdb_id,
                'protein_class': protein_class,
                'protein_description': proteins[pdb_id],
                **metrics
            })
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    return results_df

class Analysis:
    @staticmethod
    def evaluate_predictions():
        return evaluate_predictions()