#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ConSBind Tutorial - Programmatic Usage
======================================
This tutorial demonstrates how to use the ConSBind package programmatically 
for predicting protein binding sites, with a focus on analyzing the HIV-1 
protease structure (PDB ID: 1hsg).
"""

# %% [markdown]
# # ConSBind Tutorial - Programmatic Usage
# 
# This notebook demonstrates how to use the ConSBind package programmatically for predicting 
# protein binding sites, with a focus on analyzing the HIV-1 protease structure (PDB ID: 1hsg).

# %% [markdown]
# ## Overview
# 
# ConSBind (Consensus Structural Binding site predictor) identifies potential binding sites 
# in protein structures using a consensus of geometric and energy-based approaches. 
# This tutorial shows how to:
# 
# 1. Load protein structures
# 2. Run binding site predictions
# 3. Access and analyze prediction results
# 4. Visualize binding sites with advanced plots
# 5. Compare different scoring methods
# 6. Analyze pocket properties
# 
# Let's get started!

# %% [markdown]
# ## Setup
# 
# First, let's import the necessary modules and set up logging.

# %%
import os
import sys
import logging
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import Image, display, HTML
from Bio.PDB import PDBParser, Selection

# Import ConSBind modules
from ConSBind.input.file_handler import convert_ent_to_pdb, find_pdb_files
from ConSBind.core.structure import ProteinStructure
from ConSBind.core.finder import ConsensusPocketFinder
from ConSBind.core.scoring import final_scoring
from ConSBind.output.output import save_predictions, save_pymol, save_chimera

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('ConSBind-Tutorial')

# %% [markdown]
# ## Download and Prepare Protein Structure
# 
# First, we'll download the HIV-1 protease with indinavir (PDB ID: 1HSG) using BioPython's PDBList module.

# %%
import Bio.PDB.PDBList
pdblist = Bio.PDB.PDBList()
pdb_id = "1hsg" 
os.makedirs("data/tutorial", exist_ok=True)

print(f"Downloading PDB ID: {pdb_id} in PDB format...")
pdblist.retrieve_pdb_file(pdb_id, file_format="pdb", pdir="data/tutorial")
    
# PDBList downloads files with .ent extension
ent_file = Path(f"data/tutorial/pdb{pdb_id}.ent")
    
print(f"Downloaded PDB file to {ent_file}")

# Convert .ent to .pdb using ConSBind's file handler
pdb_file = str(convert_ent_to_pdb(ent_file))
print(f"Converted to PDB format: {pdb_file}")

# %% [markdown]
# ## Load Protein Structure
# 
# Now that we have a PDB file, we can load it using ConSBind's ProteinStructure class.
# This class handles parsing the PDB file and preparing it for binding site prediction.

# %%
try:
    print(f"Loading protein structure from {pdb_file}...")
    protein = ProteinStructure(pdb_file)
    print("Protein structure loaded successfully!")
    
except Exception as e:
    print(f"Error loading protein structure: {str(e)}")
    sys.exit(1)

# %% [markdown]
# ## Analyze Protein Structure
# 
# Let's examine the protein structure to understand its properties.

# %%
# Basic protein information
print(f"Protein ID: {protein.pdb_id}")
print(f"Number of chains: {len(protein.model.child_list)}")

# Count residues and atoms
residues = Selection.unfold_entities(protein.structure, 'R')
atoms = Selection.unfold_entities(protein.structure, 'A')
print(f"Number of residues: {len(residues)}")
print(f"Number of atoms: {len(atoms)}")

# %% [markdown]
# ## Visualize Protein Structure
# 
# Let's create a 3D visualization of the protein structure to better understand its shape.

# %%
def plot_protein_structure(protein, title="Protein Structure"):
    """Plot protein structure in 3D"""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get coordinates of all atoms
    atoms = Selection.unfold_entities(protein.structure, 'A')
    coords = np.array([atom.get_coord() for atom in atoms])
    
    # Plot atoms as points
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=2, alpha=0.5, c='gray')
    
    # Plot chains with different colors
    chain_colors = ['blue', 'red', 'green', 'purple', 'orange', 'cyan']
    
    for i, chain in enumerate(protein.model):
        chain_atoms = Selection.unfold_entities(chain, 'A')
        chain_coords = np.array([atom.get_coord() for atom in chain_atoms])
        
        # Use modulo to handle more chains than colors
        color = chain_colors[i % len(chain_colors)]
        
        # Plot chain backbone
        ca_atoms = [atom for atom in chain_atoms if atom.get_name() == 'CA']
        if ca_atoms:
            ca_coords = np.array([atom.get_coord() for atom in ca_atoms])
            ax.plot(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2], 
                    color=color, linewidth=2, alpha=0.7, label=f'Chain {chain.id}')
    
    # Find ligand (indinavir in 1HSG)
    hetero_atoms = [atom for atom in atoms if atom.get_id()[0].strip() not in [' ', 'H']]
    if hetero_atoms:
        hetero_coords = np.array([atom.get_coord() for atom in hetero_atoms])
        ax.scatter(hetero_coords[:, 0], hetero_coords[:, 1], hetero_coords[:, 2], 
                   s=10, c='yellow', alpha=1.0, label='Ligand')
    
    # Set labels and title
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title(title)
    
    # Add legend
    ax.legend()
    
    # Set equal aspect ratio
    max_range = np.array([
        coords[:, 0].max() - coords[:, 0].min(),
        coords[:, 1].max() - coords[:, 1].min(),
        coords[:, 2].max() - coords[:, 2].min()
    ]).max() / 2.0
    
    mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
    mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
    mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    plt.tight_layout()
    return fig

# Plot the protein structure
fig = plot_protein_structure(protein, title="HIV-1 Protease with Indinavir (1HSG)")
plt.show()

# %% [markdown]
# ## Predict Binding Sites
# 
# Now let's use ConSBind to predict binding sites in the protein structure.
# We'll use both geometric and energy-based approaches and then combine them.

# %%
# Create consensus pocket finder
pocket_finder = ConsensusPocketFinder()

# 1. Find pockets using geometric approach
print("Finding geometric pockets...")
geometric_pockets = pocket_finder.find_pockets_geometric(
    protein, 
    probe_radius=1.4,
    min_size=5
)
print(f"Found {len(geometric_pockets)} geometric pockets")

# 2. Find pockets using energy-based approach
print("Finding energy-based pockets...")
energy_pockets = pocket_finder.find_pockets_energy(
    protein,
    grid_spacing=1.0
)
print(f"Found {len(energy_pockets)} energy-based pockets")

# 3. Combine results for consensus
print("Combining results...")
consensus_pockets = pocket_finder.combine_pockets(protein, geometric_pockets, energy_pockets)
print(f"Found {len(consensus_pockets)} consensus pockets")

# 4. Score pockets
print("Scoring pockets...")
protein_type = "enzyme"  # HIV-1 protease is an enzyme
scored_pockets = final_scoring(consensus_pockets, protein_type)
print(f"Final selection: {len(scored_pockets)} significant pockets")

# %% [markdown]
# ## Analyze Binding Site Predictions
# 
# Let's analyze the predicted binding sites in more detail.

# %%
# Create a DataFrame for easier analysis
pocket_data = []
for i, pocket in enumerate(scored_pockets, 1):
    pocket_dict = {
        'Pocket': i,
        'Consensus Score': pocket['consensus_score'],
        'Final Score': pocket['final_score'],
        'Size': pocket['size'],
        'X': pocket['center'][0],
        'Y': pocket['center'][1],
        'Z': pocket['center'][2]
    }
    
    # Add other scores if available
    if 'druggability' in pocket:
        pocket_dict['Druggability'] = pocket['druggability']
    if 'knowledge_score' in pocket:
        pocket_dict['Knowledge Score'] = pocket['knowledge_score']
    if 'hydrophobicity' in pocket:
        pocket_dict['Hydrophobicity'] = pocket['hydrophobicity']
    if 'electrostatics' in pocket:
        pocket_dict['Electrostatics'] = pocket['electrostatics']
        
    pocket_data.append(pocket_dict)

pockets_df = pd.DataFrame(pocket_data)
display(pockets_df)

# %% [markdown]
# ## Visualize Binding Sites
# 
# Let's visualize the predicted binding sites in the protein structure.

# %%
def plot_protein_with_pockets(protein, pockets, title="Protein with Predicted Binding Sites"):
    """Plot protein structure with predicted binding sites"""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get coordinates of all atoms
    atoms = Selection.unfold_entities(protein.structure, 'A')
    coords = np.array([atom.get_coord() for atom in atoms])
    
    # Plot atoms as points with low alpha for transparency
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=1, alpha=0.1, c='gray')
    
    # Plot chains with different colors
    chain_colors = ['blue', 'red', 'green', 'purple', 'orange', 'cyan']
    
    for i, chain in enumerate(protein.model):
        chain_atoms = Selection.unfold_entities(chain, 'A')
        
        # Plot chain backbone
        ca_atoms = [atom for atom in chain_atoms if atom.get_name() == 'CA']
        if ca_atoms:
            ca_coords = np.array([atom.get_coord() for atom in ca_atoms])
            ax.plot(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2], 
                    color=chain_colors[i % len(chain_colors)], linewidth=1, alpha=0.5)
    
    # Find ligand (indinavir in 1HSG)
    hetero_atoms = [atom for atom in atoms if atom.get_id()[0].strip() not in [' ', 'H']]
    if hetero_atoms:
        hetero_coords = np.array([atom.get_coord() for atom in hetero_atoms])
        ax.scatter(hetero_coords[:, 0], hetero_coords[:, 1], hetero_coords[:, 2], 
                   s=20, c='yellow', alpha=1.0, label='Ligand (Indinavir)')
    
    # Plot binding sites
    pocket_colors = ['limegreen', 'orange', 'magenta', 'cyan', 'red']
    for i, pocket in enumerate(pockets):
        center = pocket['center']
        score = pocket['consensus_score']
        
        # Scale size by score
        size = 100 * score
        
        # Use modulo to handle more pockets than colors
        color = pocket_colors[i % len(pocket_colors)]
        
        ax.scatter(center[0], center[1], center[2], s=size, c=color, alpha=0.7, 
                   label=f'Pocket {i+1} (Score: {score:.2f})')
        
        # Plot pocket points if available
        if 'points' in pocket and len(pocket['points']) > 0:
            points = pocket['points']
            if isinstance(points, np.ndarray) and len(points) > 0:
                # Limit number of points to avoid cluttering
                if len(points) > 20:
                    indices = np.random.choice(len(points), 20, replace=False)
                    points = points[indices]
                
                ax.scatter(points[:, 0], points[:, 1], points[:, 2], 
                           s=10, c=color, alpha=0.3)
    
    # Set labels and title
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title(title)
    
    # Add legend
    ax.legend()
    
    # Set equal aspect ratio
    max_range = np.array([
        coords[:, 0].max() - coords[:, 0].min(),
        coords[:, 1].max() - coords[:, 1].min(),
        coords[:, 2].max() - coords[:, 2].min()
    ]).max() / 2.0
    
    mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
    mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
    mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    plt.tight_layout()
    return fig

# Plot the protein with predicted binding sites
fig = plot_protein_with_pockets(protein, scored_pockets, 
                               title="HIV-1 Protease with Predicted Binding Sites")
plt.show()

# %% [markdown]
# ## Compare Predicted Sites with Known Binding Site
# 
# The HIV-1 protease (1HSG) has a known binding site where the inhibitor indinavir binds.
# Let's compare our predictions with this known binding site.

# %%
def find_ligand_binding_site(protein):
    """Find the binding site of a ligand in the protein structure"""
    # Get all atoms
    atoms = Selection.unfold_entities(protein.structure, 'A')
    
    # Find hetero atoms (ligands)
    hetero_atoms = [atom for atom in atoms if atom.get_id()[0].strip() not in [' ', 'H']]
    
    if not hetero_atoms:
        print("No ligand found in the structure")
        return None
    
    # Calculate ligand center
    ligand_coords = np.array([atom.get_coord() for atom in hetero_atoms])
    ligand_center = np.mean(ligand_coords, axis=0)
    
    # Find protein residues near the ligand (binding site)
    binding_site_residues = []
    for residue in Selection.unfold_entities(protein.structure, 'R'):
        # Skip non-amino acid residues
        if residue.get_resname() not in protein.standard_aa_names:
            continue
        
        # Check if any atom in the residue is close to the ligand
        for atom in residue:
            min_dist = min([np.linalg.norm(atom.get_coord() - ligand_atom.get_coord()) 
                            for ligand_atom in hetero_atoms])
            if min_dist < 4.0:  # 4Å cutoff for binding site
                binding_site_residues.append(residue)
                break
    
    # Create binding site dictionary
    binding_site = {
        'center': ligand_center,
        'size': len(binding_site_residues),
        'residues': [(res.get_parent().id, res.get_id()[1], res.get_resname()) 
                     for res in binding_site_residues]
    }
    
    return binding_site

# Find the actual binding site
actual_binding_site = find_ligand_binding_site(protein)

if actual_binding_site:
    print(f"Actual binding site:")
    print(f"Center: {actual_binding_site['center']}")
    print(f"Size: {actual_binding_site['size']} residues")
    print("\nBinding site residues:")
    for chain, resid, resname in actual_binding_site['residues']:
        print(f"  {chain}:{resname}{resid}")

# %% [markdown]
# ## Compare Predictions with Actual Binding Site

# %%
def calculate_distance(point1, point2):
    """Calculate Euclidean distance between two points"""
    return np.linalg.norm(np.array(point1) - np.array(point2))

# Calculate distances between predicted pockets and actual binding site
if actual_binding_site:
    actual_center = actual_binding_site['center']
    
    distances = []
    for i, pocket in enumerate(scored_pockets, 1):
        dist = calculate_distance(pocket['center'], actual_center)
        distances.append({
            'Pocket': i,
            'Distance to Actual Site (Å)': dist,
            'Consensus Score': pocket['consensus_score'],
            'Final Score': pocket['final_score']
        })
    
    # Create DataFrame and sort by distance
    distances_df = pd.DataFrame(distances)
    distances_df = distances_df.sort_values('Distance to Actual Site (Å)')
    display(distances_df)
    
    # Plot distance vs score
    plt.figure(figsize=(10, 6))
    plt.scatter(distances_df['Distance to Actual Site (Å)'], 
                distances_df['Consensus Score'], 
                s=100, alpha=0.7)
    
    # Add labels to points
    for i, row in distances_df.iterrows():
        plt.annotate(f"Pocket {int(row['Pocket'])}", 
                     (row['Distance to Actual Site (Å)'], row['Consensus Score']),
                     xytext=(5, 5), textcoords='offset points')
    
    plt.xlabel('Distance to Actual Binding Site (Å)')
    plt.ylabel('Consensus Score')
    plt.title('Pocket Score vs Distance to Actual Binding Site')
    plt.grid(True, alpha=0.3)
    plt.show()

# %% [markdown]
# ## Save Prediction Results
# 
# Now let's save the prediction results to files for further analysis and visualization.

# %%
# Create output directory
output_dir = "results/tutorial"
os.makedirs(output_dir, exist_ok=True)
output_prefix = os.path.join(output_dir, protein.pdb_id)

# 1. Save predictions to text and PDB files
output_file, output_pdb = save_predictions(scored_pockets, protein, output_prefix)
print(f"Predictions saved to {output_file}")
print(f"Modified PDB saved to {output_pdb}")

# 2. Generate PyMOL visualization script
pymol_script = save_pymol(scored_pockets, protein, output_prefix, output_pdb)
print(f"PyMOL script saved to {pymol_script}")

# 3. Generate Chimera visualization script
chimera_script = save_chimera(scored_pockets, protein, output_prefix, output_pdb)
print(f"Chimera script saved to {chimera_script}")

# %% [markdown]
# ## Visualize Binding Site Properties
# 
# Let's analyze the properties of the predicted binding sites in more detail.

# %%
# Extract residues for each pocket
pocket_residues = []
for i, pocket in enumerate(scored_pockets, 1):
    residues = protein.get_pocket_residues(pocket)
    pocket_residues.append({
        'Pocket': i,
        'Residues': residues,
        'Score': pocket['consensus_score']
    })

# Calculate amino acid composition for each pocket
aa_counts = []
for pocket in pocket_residues:
    # Count amino acids by type
    aa_count = {}
    for _, _, resname in pocket['Residues']:
        if resname in aa_count:
            aa_count[resname] += 1
        else:
            aa_count[resname] = 1
    
    # Add pocket info
    aa_count['Pocket'] = pocket['Pocket']
    aa_count['Score'] = pocket['Score']
    aa_counts.append(aa_count)

# Convert to DataFrame
aa_df = pd.DataFrame(aa_counts)
aa_df = aa_df.fillna(0)  # Replace NaN with 0

# Plot amino acid composition
plt.figure(figsize=(14, 8))
for i, pocket in enumerate(aa_counts):
    pocket_id = pocket['Pocket']
    
    # Extract amino acid counts (exclude Pocket and Score columns)
    aa_names = [col for col in aa_df.columns if col not in ['Pocket', 'Score']]
    aa_values = [pocket.get(aa, 0) for aa in aa_names]
    
    # Plot as bar chart
    plt.subplot(1, len(aa_counts), i+1)
    plt.bar(aa_names, aa_values)
    plt.title(f'Pocket {pocket_id} (Score: {pocket["Score"]:.2f})')
    plt.xticks(rotation=90)
    plt.tight_layout()

plt.suptitle('Amino Acid Composition of Predicted Binding Sites', y=1.05)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Programmatic Access to Pocket Data
# 
# Let's demonstrate how to access and work with the pocket data programmatically.

# %%
# Example: Find the highest scoring pocket
highest_scoring = max(scored_pockets, key=lambda x: x['consensus_score'])
print(f"Highest scoring pocket:")
print(f"Consensus Score: {highest_scoring['consensus_score']:.2f}")
print(f"Final Score: {highest_scoring['final_score']:.2f}")
print(f"Size: {highest_scoring['size']}")
print(f"Center: {highest_scoring['center']}")

# Example: Calculate average properties of all pockets
avg_consensus = np.mean([p['consensus_score'] for p in scored_pockets])
avg_size = np.mean([p['size'] for p in scored_pockets])
print(f"\nAverage properties:")
print(f"Average Consensus Score: {avg_consensus:.2f}")
print(f"Average Size: {avg_size:.2f}")

# Example: Find pockets with specific properties
large_pockets = [p for p in scored_pockets if p['size'] > 100]
print(f"\nNumber of large pockets (size > 100): {len(large_pockets)}")

high_scoring = [p for p in scored_pockets if p['consensus_score'] > 0.7]
print(f"Number of high-scoring pockets (score > 0.7): {len(high_scoring)}")

# %% [markdown]
# ## Additional Example: Comparing Apo and Holo Forms of a Protein
# 
# A common analysis in structural biology is comparing the binding sites of a protein in its apo (unbound) and holo (ligand-bound) forms. This comparison can reveal conformational changes that occur upon ligand binding.
# 
# For this example, we'll use adenylate kinase (AdK), which undergoes significant conformational changes upon binding:
# - Apo form: 4AKE (open conformation)
# - Holo form: 1AKE (closed conformation with inhibitor)

# %%
# Download and prepare both structures
print("Downloading and preparing apo and holo structures...")

# Download apo structure (4AKE)
apo_pdb_id = "4ake"
pdblist.retrieve_pdb_file(apo_pdb_id, file_format="pdb", pdir="data/tutorial")
apo_ent_file = Path(f"data/tutorial/pdb{apo_pdb_id}.ent")
apo_pdb_file = str(convert_ent_to_pdb(apo_ent_file))
print(f"Apo structure (4AKE) prepared: {apo_pdb_file}")

# Download holo structure (1AKE)
holo_pdb_id = "1ake"
pdblist.retrieve_pdb_file(holo_pdb_id, file_format="pdb", pdir="data/tutorial")
holo_ent_file = Path(f"data/tutorial/pdb{holo_pdb_id}.ent")
holo_pdb_file = str(convert_ent_to_pdb(holo_ent_file))
print(f"Holo structure (1AKE) prepared: {holo_pdb_file}")

# %%
# Load both structures
apo_protein = ProteinStructure(apo_pdb_file)
holo_protein = ProteinStructure(holo_pdb_file)

print("Structures loaded successfully!")
print(f"Apo structure (4AKE): {len(Selection.unfold_entities(apo_protein.structure, 'R'))} residues")
print(f"Holo structure (1AKE): {len(Selection.unfold_entities(holo_protein.structure, 'R'))} residues")

# %%
# Predict binding sites for both structures
pocket_finder = ConsensusPocketFinder()

# Process apo structure
print("\nProcessing apo structure (4AKE)...")
apo_geometric = pocket_finder.find_pockets_geometric(apo_protein)
apo_energy = pocket_finder.find_pockets_energy(apo_protein)
apo_consensus = pocket_finder.combine_pockets(apo_protein, apo_geometric, apo_energy)
apo_pockets = final_scoring(apo_consensus, "enzyme")  # AdK is an enzyme
print(f"Found {len(apo_pockets)} significant pockets in apo structure")

# Process holo structure
print("\nProcessing holo structure (1AKE)...")
holo_geometric = pocket_finder.find_pockets_geometric(holo_protein)
holo_energy = pocket_finder.find_pockets_energy(holo_protein)
holo_consensus = pocket_finder.combine_pockets(holo_protein, holo_geometric, holo_energy)
holo_pockets = final_scoring(holo_consensus, "enzyme")
print(f"Found {len(holo_pockets)} significant pockets in holo structure")

# %%
# Find the actual binding site in the holo structure
holo_binding_site = find_ligand_binding_site(holo_protein)

if holo_binding_site:
    print("\nActual binding site in holo structure (1AKE):")
    print(f"Center: {holo_binding_site['center']}")
    print(f"Size: {holo_binding_site['size']} residues")

# %%
# Visualize both structures with their predicted binding sites
fig, axes = plt.subplots(1, 2, figsize=(20, 10), subplot_kw={'projection': '3d'})

# Plot apo structure
ax = axes[0]
atoms = Selection.unfold_entities(apo_protein.structure, 'A')
coords = np.array([atom.get_coord() for atom in atoms])
ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=1, alpha=0.1, c='gray')

# Plot apo binding sites
for i, pocket in enumerate(apo_pockets):
    center = pocket['center']
    score = pocket['consensus_score']
    size = 100 * score
    color = ['limegreen', 'orange', 'magenta', 'cyan', 'red'][i % 5]
    ax.scatter(center[0], center[1], center[2], s=size, c=color, alpha=0.7, 
               label=f'Pocket {i+1} (Score: {score:.2f})')

ax.set_title('Apo Structure (4AKE) with Predicted Binding Sites')
ax.legend()

# Plot holo structure
ax = axes[1]
atoms = Selection.unfold_entities(holo_protein.structure, 'A')
coords = np.array([atom.get_coord() for atom in atoms])
ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=1, alpha=0.1, c='gray')

# Plot holo binding sites
for i, pocket in enumerate(holo_pockets):
    center = pocket['center']
    score = pocket['consensus_score']
    size = 100 * score
    color = ['limegreen', 'orange', 'magenta', 'cyan', 'red'][i % 5]
    ax.scatter(center[0], center[1], center[2], s=size, c=color, alpha=0.7, 
               label=f'Pocket {i+1} (Score: {score:.2f})')

# Plot actual binding site if available
if holo_binding_site:
    center = holo_binding_site['center']
    ax.scatter(center[0], center[1], center[2], s=200, c='yellow', alpha=0.7, 
               label='Actual Binding Site')

ax.set_title('Holo Structure (1AKE) with Predicted Binding Sites')
ax.legend()

plt.tight_layout()
plt.show()

# %%
# Compare binding site residues between apo and holo forms
def get_residue_ids(protein, pocket):
    """Get standardized residue IDs for comparison"""
    residues = protein.get_pocket_residues(pocket)
    return set([f"{chain}:{resname}{resid}" for chain, resid, resname in residues])

# Create DataFrames for comparison
apo_data = []
for i, pocket in enumerate(apo_pockets, 1):
    residue_ids = get_residue_ids(apo_protein, pocket)
    apo_data.append({
        'Pocket': f"Apo-{i}",
        'Score': pocket['consensus_score'],
        'Size': pocket['size'],
        'Residue_IDs': residue_ids
    })

holo_data = []
for i, pocket in enumerate(holo_pockets, 1):
    residue_ids = get_residue_ids(holo_protein, pocket)
    holo_data.append({
        'Pocket': f"Holo-{i}",
        'Score': pocket['consensus_score'],
        'Size': pocket['size'],
        'Residue_IDs': residue_ids
    })

# Calculate overlap between pockets
overlap_data = []
for apo_entry in apo_data:
    for holo_entry in holo_data:
        apo_residues = apo_entry['Residue_IDs']
        holo_residues = holo_entry['Residue_IDs']
        
        # Calculate overlap
        common_residues = apo_residues.intersection(holo_residues)
        overlap_percent = len(common_residues) / len(apo_residues) * 100 if apo_residues else 0
        
        overlap_data.append({
            'Apo_Pocket': apo_entry['Pocket'],
            'Holo_Pocket': holo_entry['Pocket'],
            'Apo_Score': apo_entry['Score'],
            'Holo_Score': holo_entry['Score'],
            'Common_Residues': len(common_residues),
            'Overlap_Percent': overlap_percent
        })

# Convert to DataFrame and display
overlap_df = pd.DataFrame(overlap_data)
overlap_df = overlap_df.sort_values('Overlap_Percent', ascending=False)
display(overlap_df.head(10))

# %%
# Visualize the overlap between apo and holo pockets
plt.figure(figsize=(12, 8))
sns.heatmap(
    overlap_df.pivot(index='Apo_Pocket', columns='Holo_Pocket', values='Overlap_Percent'),
    annot=True, 
    cmap='YlGnBu',
    fmt='.1f'
)
plt.title('Overlap Between Apo and Holo Binding Sites (%)')
plt.tight_layout()
plt.show()

# %%
# Analyze conformational changes between apo and holo forms
# Find the best matching pockets (highest overlap)
best_match = overlap_df.iloc[0]
apo_pocket_idx = int(best_match['Apo_Pocket'].split('-')[1]) - 1
holo_pocket_idx = int(best_match['Holo_Pocket'].split('-')[1]) - 1

apo_pocket = apo_pockets[apo_pocket_idx]
holo_pocket = holo_pockets[holo_pocket_idx]

print(f"Best matching pockets: {best_match['Apo_Pocket']} and {best_match['Holo_Pocket']}")
print(f"Overlap: {best_match['Overlap_Percent']:.1f}% ({best_match['Common_Residues']} common residues)")

# Calculate distance between pocket centers
distance = calculate_distance(apo_pocket['center'], holo_pocket['center'])
print(f"Distance between pocket centers: {distance:.2f} Å")

# Calculate volume change
volume_change = holo_pocket['size'] - apo_pocket['size']
volume_percent = volume_change / apo_pocket['size'] * 100
print(f"Volume change: {volume_change:.1f} units ({volume_percent:.1f}%)")

# %% [markdown]
# ## Key Insights from Apo vs. Holo Comparison
# 
# The comparison between the apo (4AKE) and holo (1AKE) forms of adenylate kinase reveals:
# 
# 1. **Conformational Changes**: The protein undergoes significant conformational changes upon ligand binding, with the binding pocket becoming more defined in the holo form.
# 
# 2. **Binding Site Conservation**: Despite these changes, many of the same residues are involved in forming the binding site in both conformations.
# 
# 3. **Pocket Volume**: The binding pocket often changes in volume, typically becoming smaller and more tightly packed around the ligand in the holo form.
# 
# 4. **Prediction Accuracy**: ConSBind can identify potential binding sites in both conformations, with higher scores typically assigned to the more defined pocket in the holo form.
# 
# This type of analysis is valuable for understanding protein flexibility and how it relates to function, as well as for drug design where targeting specific conformational states may be desirable.

# %% [markdown]
# ## Conclusion
# 
# In this tutorial, we've demonstrated how to use ConSBind programmatically to:
# 
# 1. Load a protein structure (HIV-1 protease, PDB ID: 1hsg)
# 2. Predict binding sites using geometric and energy-based approaches
# 3. Combine and score the predictions
# 4. Analyze and visualize the results
# 5. Compare predictions with the known binding site
# 6. Save the results for further analysis
# 7. Compare apo and holo forms of a protein to analyze conformational changes
# 
# ConSBind provides a powerful and flexible framework for binding site prediction that can be
# integrated into larger workflows for drug discovery, protein function analysis, and more.