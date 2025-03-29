#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein Structure Module
========================
This module handles the loading and analysis of protein structures.
"""

import os
import numpy as np
import logging
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from Bio.PDB import PDBParser, Selection, PDBIO, NeighborSearch
from Bio.PDB.DSSP import DSSP

# Use the same logger as main to prevent duplicate messages
logger = logging.getLogger('ConSBind')

class ProteinStructure:
    """Class to handle protein structure analysis"""
    
    def __init__(self, pdb_file):
        """Initialize with a PDB file"""
        self.pdb_file = pdb_file
        self.pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
        
        # Parse PDB file
        parser = PDBParser(QUIET=True)
        try:
            self.structure = parser.get_structure(self.pdb_id, pdb_file)
            self.model = self.structure[0]  # Use the first model
        except Exception as e:
            logger.error(f"Failed to parse PDB file: {e}")
            raise
            
        # Calculate structure properties
        self.calculate_surface_properties()
        
    def calculate_surface_properties(self):
        """Calculate surface properties using DSSP"""
        try:
            # Run DSSP to get accessible surface area
            dssp = DSSP(self.model, self.pdb_file, dssp='mkdssp')
            self.dssp_data = dssp
        except Exception as e:
            logger.warning(f"DSSP calculation failed: {e}")
            self.dssp_data = None
    
    def get_surface_atoms(self, rel_asa_threshold=0.2):
        """Get atoms on the protein surface based on relative accessible surface area"""
        surface_atoms = []
        
        if self.dssp_data is None:
            # If DSSP failed, use distance-based approach
            all_atoms = Selection.unfold_entities(self.structure, 'A')
            
            # Create KDTree for efficient neighbor search
            coords = np.array([atom.get_coord() for atom in all_atoms])
            kdtree = KDTree(coords)
            
            # Identify surface atoms as those with fewer neighbors
            for i, atom in enumerate(all_atoms):
                # Count neighbors within 8Å
                neighbors = kdtree.query_ball_point(atom.get_coord(), 8.0)
                if len(neighbors) < 15:  # Threshold for surface atoms
                    surface_atoms.append(atom)
        else:
            # Use DSSP data to identify surface residues
            surface_residues = []
            for key in self.dssp_data.keys():
                data = self.dssp_data[key]
                rel_asa = data[3]  # Relative accessible surface area
                if rel_asa > rel_asa_threshold:
                    chain_id, residue_id = key[0], key[1]
                    residue = self.model[chain_id][residue_id]
                    surface_residues.append(residue)
            
            # Get atoms from surface residues
            for residue in surface_residues:
                for atom in residue:
                    if atom.get_name() not in ['H', 'HA']:  # Skip hydrogen atoms
                        surface_atoms.append(atom)
        
        logger.info(f"Identified {len(surface_atoms)} surface atoms")
        return surface_atoms
    
    def get_cavities(self, probe_radius=1.4, grid_spacing=1.0, min_cavity_size=5, detect_filled = True):
        """
        Find cavities using a grid-based approach, with option to detect filled cavities
        """
        # Get protein atoms
        atoms = Selection.unfold_entities(self.model, 'A')
        coords = np.array([atom.get_coord() for atom in atoms])

        # Identify possible hetero atoms and exclude them for cavity detection 
        hetero_coords = None 
        if detect_filled:
            hetero_atoms = [atom for atom in atoms if atom.get_id()[0].strip() not in [' ', 'H']]
            if hetero_atoms:
                hetero_coords = np.array([atom.get_coord() for atom in hetero_atoms])
        
        # Define grid around the protein
        padding = 10.0  # Å
        min_coords = np.min(coords, axis=0) - padding
        max_coords = np.max(coords, axis=0) + padding
        
        # Generate grid points
        x = np.arange(min_coords[0], max_coords[0], grid_spacing)
        y = np.arange(min_coords[1], max_coords[1], grid_spacing)
        z = np.arange(min_coords[2], max_coords[2], grid_spacing)
        
        logger.info(f"Created grid with dimensions: {len(x)}x{len(y)}x{len(z)}")
        
        # Create KDTree for efficient distance calculations
        kdtree = KDTree(coords)
        
        # Identify cavity points
        cavity_points = []
        
        # Sample grid points
        sample_size = min(10000, len(x) * len(y) * len(z))
        indices = np.random.choice(len(x) * len(y) * len(z), size=sample_size, replace=False)
        
        for idx in indices:
            # Convert linear index to 3D coordinates
            i = idx // (len(y) * len(z))
            j = (idx % (len(y) * len(z))) // len(z)
            k = idx % len(z)
            
            if i >= len(x) or j >= len(y) or k >= len(z):
                continue
                
            point = np.array([x[i], y[j], z[k]])
            
            # Find distance to nearest atom
            dist, _ = kdtree.query(point)

            # Check if point is inside protein but not too close to atoms
            cavity_condition = probe_radius < dist < 4.0
        
            # If detect_filled is True, also consider points near hetero atoms as potential cavities
            if detect_filled and hetero_coords is not None and not cavity_condition:
                # Calculate distance to nearest hetero atom
                hetero_dists = np.linalg.norm(hetero_coords - point, axis=1)
                min_hetero_dist = np.min(hetero_dists) if len(hetero_dists) > 0 else float('inf')
                
                # Points near hetero atoms could be binding sites
                if min_hetero_dist < 4.0:
                    cavity_condition = True
            
            if cavity_condition:
                # Check if the point is enclosed by protein atoms
                directions = np.array([
                    [1, 0, 0], [-1, 0, 0], 
                    [0, 1, 0], [0, -1, 0],
                    [0, 0, 1], [0, 0, -1]
                ])
                
                enclosed = True
                for direction in directions:
                    ray_point = point
                    max_dist = 10.0  # Maximum ray distance
                    hit_protein = False
                    
                    # Cast ray
                    for t in np.arange(1.0, max_dist, 1.0):
                        ray_point = point + t * direction
                        # Check if ray hits protein
                        ray_dist, _ = kdtree.query(ray_point)
                        if ray_dist < probe_radius:
                            hit_protein = True
                            break
                    
                    if not hit_protein:
                        enclosed = False
                        break
                
                if enclosed:
                    cavity_points.append(point)
            
            # Check if point is inside protein but not too close to atoms
            if probe_radius < dist < 4.0:  # Inside a pocket but not inside an atom
                # Check if the point is enclosed by protein atoms
                directions = np.array([
                    [1, 0, 0], [-1, 0, 0], 
                    [0, 1, 0], [0, -1, 0],
                    [0, 0, 1], [0, 0, -1]
                ])
                
                enclosed = True
                for direction in directions:
                    ray_point = point
                    max_dist = 10.0  # Maximum ray distance
                    hit_protein = False
                    
                    # Cast ray
                    for t in np.arange(1.0, max_dist, 1.0):
                        ray_point = point + t * direction
                        # Check if ray hits protein
                        ray_dist, _ = kdtree.query(ray_point)
                        if ray_dist < probe_radius:
                            hit_protein = True
                            break
                    
                    if not hit_protein:
                        enclosed = False
                        break
                
                if enclosed:
                    cavity_points.append(point)
        
        logger.info(f"Found {len(cavity_points)} potential cavity points")
        
        # Cluster cavity points
        if len(cavity_points) > 1:
            cavity_points = np.array(cavity_points)
            distances = pdist(cavity_points)
            linkage_matrix = linkage(distances, method='single')
            clusters = fcluster(linkage_matrix, t=3.0, criterion='distance')
            
            # Filter clusters by size
            unique_clusters, counts = np.unique(clusters, return_counts=True)
            valid_clusters = unique_clusters[counts >= min_cavity_size]
            
            # Calculate cluster centers
            cavity_centers = []
            for cluster_id in valid_clusters:
                cluster_mask = clusters == cluster_id
                cluster_points = cavity_points[cluster_mask]
                center = np.mean(cluster_points, axis=0)
                cavity_centers.append({
                    'center': center,
                    'size': len(cluster_points),
                    'points': cluster_points
                })
            
            logger.info(f"Found {len(cavity_centers)} cavities after clustering")
            return cavity_centers
        else:
            return []

    def calculate_hydrophobicity(self, center, radius=8.0):
        """Calculate average hydrophobicity around a center point"""
        # Kyte & Doolittle hydrophobicity scale
        hydrophobicity = {
            'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5, 'MET': 1.9, 'ALA': 1.8,
            'GLY': -0.4, 'THR': -0.7, 'SER': -0.8, 'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6,
            'HIS': -3.2, 'GLU': -3.5, 'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
        }
        
        # Find residues within radius of center
        atoms = Selection.unfold_entities(self.model, 'A')
        ns = NeighborSearch(atoms)
        nearby_atoms = ns.search(center, radius, 'R')  # Search for residues
        
        if not nearby_atoms:
            return 0.0
        
        # Calculate average hydrophobicity
        total = 0
        count = 0
        for residue in nearby_atoms:
            if residue.get_resname() in hydrophobicity:
                total += hydrophobicity[residue.get_resname()]
                count += 1
        
        return total / max(1, count)

    def calculate_electrostatics(self, center, radius=8.0):
        """Calculate an electrostatic score around the center point"""
        # Simple charge assignment
        charges = {
            'ARG': 1.0, 'LYS': 1.0, 'HIS': 0.5,  # Positive
            'ASP': -1.0, 'GLU': -1.0,  # Negative
            'SER': 0.1, 'THR': 0.1, 'ASN': 0.1, 'GLN': 0.1, 'TYR': 0.1  # Polar
        }
        
        # Find residues within radius of center
        atoms = Selection.unfold_entities(self.model, 'A')
        ns = NeighborSearch(atoms)
        nearby_atoms = ns.search(center, radius, 'R')  # Search for residues
        
        # Sum charges with distance weighting
        total_charge = 0.0
        for residue in nearby_atoms:
            if residue.get_resname() in charges:
                # Get closest atom in residue to center
                min_dist = float('inf')
                for atom in residue:
                    dist = np.linalg.norm(atom.get_coord() - center)
                    min_dist = min(min_dist, dist)
                
                # Apply distance-weighted charge
                charge = charges[residue.get_resname()]
                weight = 1.0 / max(1.0, min_dist)  # Simple distance weighting
                total_charge += charge * weight
        
        return total_charge
    
    def get_pocket_residues(self, pocket, radius=8.0):
        """Get residues within a certain radius of a pocket center"""
        atoms = Selection.unfold_entities(self.model, 'A')
        ns = NeighborSearch(atoms)
        nearby_atoms = ns.search(pocket['center'], radius, 'R')  # Search for residues
        
        # Create a list of unique residues
        residues = set()
        for residue in nearby_atoms:
            # Only consider standard amino acids
            if residue.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                                    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
                                    'TYR', 'VAL']:
                chain_id = residue.get_parent().id
                res_id = residue.id[1]
                res_name = residue.get_resname()
                residues.add((chain_id, res_id, res_name))
        
        return sorted(residues)
