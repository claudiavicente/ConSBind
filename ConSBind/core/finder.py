#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Consensus Pocket Finder Module
==============================
This module implements different pocket detection methods and
combines their results for consensus-based binding site prediction.
"""

import numpy as np
import logging
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import KDTree
from Bio.PDB import Selection

# Get logger but prevent duplicate messages
logger = logging.getLogger('ConSBind')

class ConsensusPocketFinder:
    """Class that implements different pocket detection methods and combines their results"""
    
    def __init__(self):
        """Initialize the pocket finder with default parameters"""
        pass
    
    def find_pockets_geometric(self, protein, probe_radius=1.4, min_size=5):
        """Find pockets using geometric approach"""
        # Cavity detection 
        cavities = protein.get_cavities(probe_radius=probe_radius, min_cavity_size=min_size)

        # If no cavities found, try more aggressive parameters 
        if not cavities:
            logger.info("No cavities found with default parameters, trying alternatives...")

            # Try larger prove radius for larger cavities 
            cavities = protein.get_cavities(probe_radius=1.8, min_cavity_size=3)

            # If still no cavities, try surface-based approach:
            if not cavities:
                cavities = self._find_surface_pockets(protein)
        
        return cavities
    
    def _find_surface_pockets(self, protein):
        """ Alternative method to find potential binding sites in surface contours"""
        surface_atoms = protein.get_surface_atoms(rel_asa_threshold=0.15)

        # Use clustering to identify potential pocket regions 
        if len(surface_atoms) > 5:
            coords = np.array([atom.get_coord() for atom in surface_atoms])

            # Use DBSCAN for clustering 
            from sklearn.cluster import DBSCAN # type: ignore
            clustering = DBSCAN(eps=3.5, min_samples=5).fit(coords)

            labels = clustering.labels_
            unique_labels = set(labels)

            cavities = []
            for label in unique_labels:
                if label == -1:     # Skip noise 
                    continue 

                cluster_points = coords[labels == label]
                center = np.mean(cluster_points, axis=0)

                # Check if this is a pocket-like feature (concave region on a protein surface)
                if self._is_concave(protein, center):
                    cavities.append({
                        'center': center, 
                        'size': len(cluster_points),
                        'points': cluster_points
                    })
            
            return cavities
        
        return []
    
    def _is_concave(self, protein, point, num_directions=20):
        """Check if a point is in a concave region by ray-casting"""
        # Generate random directions 
        directions = np.random.randn(num_directions, 3)
        directions = directions / np.linalg.norm(directions, axis=1)[:, np.newaxis]

        # Cast rays and check for protein hits 
        atoms = Selection.unfold_entities(protein.model, 'A')
        coords = np.array([atom.get_coord() for atom in atoms])
        kdtree = KDTree(coords)

        hit_count = 0
        for direction in directions:
            hit = False
            for dist in range(1, 10):
                test_point = point + dist * direction
                nearest_dist, _ = kdtree.query(test_point)
                if nearest_dist < 2.0:  # Hit protein
                    hit = True 
                    break
            if hit:
                hit_count += 1
        
        # If most rays hit protein, it's likely concave 
        return hit_count >= num_directions * 0.7 
    
    def find_pockets_energy(self, protein, grid_spacing=1.0):
        """
        Find pockets using energy-based approach
        This is a simplified version focusing on hydrophobicity and charge
        """
        # Get protein surface
        surface_atoms = protein.get_surface_atoms()
        if not surface_atoms:
            logger.warning("No surface atoms found")
            return []
        
        # Create initial grid around the protein
        coords = np.array([atom.get_coord() for atom in surface_atoms])
        min_coords = np.min(coords, axis=0) - 5.0
        max_coords = np.max(coords, axis=0) + 5.0
        
        # Create a sparse grid (for performance)
        x = np.arange(min_coords[0], max_coords[0], grid_spacing * 2)
        y = np.arange(min_coords[1], max_coords[1], grid_spacing * 2)
        z = np.arange(min_coords[2], max_coords[2], grid_spacing * 2)
        
        # Calculate energy score for each grid point
        energy_points = []
        
        # Sample grid points
        logger.info("Calculating energy scores for grid points...")
        sample_size = min(5000, len(x) * len(y) * len(z))
        indices = np.random.choice(len(x) * len(y) * len(z), size=sample_size, replace=False)
        
        for idx in indices:
            # Convert linear index to 3D coordinates
            i = idx // (len(y) * len(z))
            j = (idx % (len(y) * len(z))) // len(z)
            k = idx % len(z)
            
            if i >= len(x) or j >= len(y) or k >= len(z):
                continue
                
            point = np.array([x[i], y[j], z[k]])
            
            # Calculate distance to nearest atom
            min_dist = float('inf')
            for atom in surface_atoms:
                dist = np.linalg.norm(atom.get_coord() - point)
                min_dist = min(min_dist, dist)
            
            # Skip points too far from or too close to the protein
            if min_dist < 1.0 or min_dist > 5.0:
                continue
            
            # Calculate scores for this point
            hydrophobicity = protein.calculate_hydrophobicity(point)
            electrostatics = protein.calculate_electrostatics(point)
            
            # Combined score - higher is better
            # Favorable pockets are usually hydrophobic with moderate charge
            energy_score = (hydrophobicity * 2.0 + 
                           abs(electrostatics) * 0.5)
            
            if energy_score > 3.0:  # Threshold for potential binding sites
                energy_points.append({
                    'center': point,
                    'score': energy_score,
                    'hydrophobicity': hydrophobicity,
                    'electrostatics': electrostatics
                })
        
        logger.info(f"Found {len(energy_points)} high-energy points")
        
        # Cluster energy points
        if len(energy_points) > 1:
            points = np.array([p['center'] for p in energy_points])
            distances = pdist(points)
            linkage_matrix = linkage(distances, method='average')
            clusters = fcluster(linkage_matrix, t=3.5, criterion='distance')
            
            # Calculate cluster centers and scores
            unique_clusters = np.unique(clusters)
            energy_clusters = []
            
            for cluster_id in unique_clusters:
                cluster_mask = clusters == cluster_id
                cluster_indices = np.where(cluster_mask)[0]
                
                if len(cluster_indices) < 3:  # Skip very small clusters
                    continue
                
                cluster_points = points[cluster_mask]
                center = np.mean(cluster_points, axis=0)
                
                # Calculate average score
                score = np.mean([energy_points[i]['score'] for i in cluster_indices])
                
                energy_clusters.append({
                    'center': center,
                    'points': cluster_points,
                    'size': len(cluster_points),
                    'score': score
                })
            
            logger.info(f"Found {len(energy_clusters)} energy-based pockets after clustering")
            return energy_clusters
        else:
            return []
    
    def combine_pockets(self, protein, geometric_pockets, energy_pockets, distance_threshold=5.0):
        """Combine pockets from different methods"""
        if not geometric_pockets and not energy_pockets:
            return []
        
        # Start with all pockets
        all_pockets = []
        
        # Add geometric pockets
        for pocket in geometric_pockets:
            # Calculate additional scores 
            druggability = self.calculate_druggability_score(protein, pocket)
            knowledge_score = self.evaluate_with_knowledge_base(protein, pocket)

            # Combined score (weighted)
            total_score = (
                pocket['size'] * 0.2 +  # Geometric size
                druggability * 1.0 +    # Druggability
                knowledge_score * 1.2   # Knowledge-based score 
            )

            all_pockets.append({
                'center': pocket['center'],
                'size': pocket['size'],
                'method': 'geometric',
                'methods': ['geometric'],  # Track all methods that detected this pocket
                'druggability': druggability, 
                'knowledge_score': knowledge_score, 
                'score': total_score,  
                'consensus_score': 1,  # Start with base score
                'points': pocket.get('points', [])
            })
        
        # Add energy pockets
        for pocket in energy_pockets:
            all_pockets.append({
                'center': pocket['center'],
                'size': pocket['size'],
                'method': 'energy',
                'methods': ['energy'],  # Track all methods that detected this pocket
                'score': pocket['score'],
                'consensus_score': 1,  # Start with base score
                'points': pocket.get('points', [])
            })
        
        # Increase consensus score for pockets that are close to each other
        for i in range(len(all_pockets)):
            for j in range(i+1, len(all_pockets)):
                dist = np.linalg.norm(all_pockets[i]['center'] - all_pockets[j]['center'])
                if dist < distance_threshold:
                    all_pockets[i]['consensus_score'] += 2  # Increased weight for consensus
                    all_pockets[j]['consensus_score'] += 2
                    
                    # Add methods from j to i's methods list
                    all_pockets[i]['methods'].extend([m for m in all_pockets[j]['methods'] if m not in all_pockets[i]['methods']])
                    # Add methods from i to j's methods list
                    all_pockets[j]['methods'].extend([m for m in all_pockets[i]['methods'] if m not in all_pockets[j]['methods']])

        # Re-rank by combined score, with consensus as primary criterion
        for pocket in all_pockets:
            pocket['final_score'] = (
                pocket['consensus_score'] * 3.0 +   # Primary criterion, so higher weight
                pocket.get('score', 0) * 0.5        # Secondary criterion with reduced weight
            )
        
        # Sort by final score
        all_pockets.sort(key=lambda x: x['final_score'], reverse=True)
        
        # Filter overlapping pockets
        filtered_pockets = []
        for pocket in all_pockets:
            # Check if this pocket is too close to any already selected pocket
            overlap = False
            for selected in filtered_pockets:
                dist = np.linalg.norm(pocket['center'] - selected['center'])
                if dist < distance_threshold:
                    overlap = True
                    # If this pocket was found by different methods than the selected one,
                    # add those methods to the selected pocket's methods list
                    for method in pocket['methods']:
                        if method not in selected['methods']:
                            selected['methods'].append(method)
                    break
            
            if not overlap:
                filtered_pockets.append(pocket)
        
        logger.info(f"Combined into {len(filtered_pockets)} consensus pockets")
        return filtered_pockets
    
    # Knowledge-based filtering 
    def evaluate_with_knowledge_base(self, protein, pocket):
        """Score pocket based on binding site knowledge base"""
        # Get residues in the pocket 
        pocket_residues = protein.get_pocket_residues(pocket)
        residue_types = [res[2] for res in pocket_residues]     # Extract residue names

        # Features that caracterize binding sites 
        score = 0.0

        # 1. Binding sites typically have a mix of hydrophobic and polar residues 
        hydrophobic = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']
        polar = ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN', 'HIS']
        charged = ['LYS', 'ARG', 'ASP', 'GLU']

        hydrophobic_count = sum(1 for res in residue_types if res in hydrophobic)
        polar_count = sum(1 for res in residue_types if res in polar)
        charged_count = sum(1 for res in residue_types if res in charged)

        # Binding sites typically have a good mix of residue types 
        if hydrophobic_count > 0 and (polar_count > 0 or charged_count > 0):
            score += 1.0

        # 2. Presence of key catalytic residues (common in enzyme active sites)
        catalytic_pairs = [
            {'HIS', 'ASP'}, {'SER', 'HIS'}, {'CYS', 'HIS'}, 
            {'LYS', 'ASP'}, {'ARG', 'ASP'}, {'ARG', 'GLU'}
        ]

        residue_set = set(residue_types)
        for pair in catalytic_pairs:
            if pair.issubset(residue_set):
                score += 0.5
        
        # 3. Presence of aromatic residues (common in binding sites)
        aromatic_count = sum(1 for res in residue_types if res in ['PHE', 'TYR', 'TRP', 'HIS'])
        if aromatic_count >= 2:
            score += 0.5

        # 4. Specific binding site characteristics
    
        # Heme binding sites often have HIS, MET, CYS
        if 'HIS' in residue_set and ('MET' in residue_set or 'CYS' in residue_set):
            score += 1.0
        
        # Nucleotide binding sites often have a glycine-rich loop and charged residues
        if residue_types.count('GLY') >= 3 and ('LYS' in residue_set or 'ARG' in residue_set):
            score += 1.0
        
        # Metal binding sites often have specific coordinating residues
        metal_binding = ['HIS', 'CYS', 'ASP', 'GLU']
        if sum(1 for res in residue_types if res in metal_binding) >= 3:
            score += 1.0
        
        # Substrate binding in enzymes - look for catalytic triad-like patterns
        if {'SER', 'HIS', 'ASP'}.issubset(residue_set) or {'CYS', 'HIS', 'ASP'}.issubset(residue_set):
            score += 1.5

        return score 

    # Druggability score 
    def calculate_druggability_score(self, protein, pocket):
        """Calculate a druggability score for the pocket"""
        # Get residues in the pocket 
        residues = protein.get_pocket_residues(pocket)

        # Calculate volume and hydrophobicity
        volume = pocket['size'] * 8.0   # Approximate
        hydrophobicity = protein.calculate_hydrophobicity(pocket['center'])

        # Mixture of hydrophobic and hydrophilic residues 
        residue_types = [res[2] for res in residues]
        hydrophobic = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']
        polar = ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN', 'HIS', 'LYS', 'ARG', 'ASP', 'GLU']

        hydrophobic_frac = sum(1 for res in residue_types if res in hydrophobic) / max(1, len(residue_types))

        # Druggability score - empirically derived. Favorable pockets have volume 200-800A, moderate hydrophobicity (0.4-0.8)
        volume_score = max(0, 1 - abs(volume - 500) / 300)
        hydrophobic_score = max(0, 1 - abs(hydrophobic_frac - 0.6) / 0.4)

        druggability = (volume_score + hydrophobic_score) / 2.0

        return druggability
    