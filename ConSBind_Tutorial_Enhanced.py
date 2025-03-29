#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# ConSBind Tutorial: Programmatic Usage (Enhanced)

This notebook demonstrates how to use the ConSBind package programmatically for predicting protein binding sites,
with enhanced visualizations and analysis.

## Overview

ConSBind (Consensus Structural Binding site predictor) identifies potential binding sites in protein structures
using a consensus of geometric and energy-based approaches. This tutorial shows how to:

1. Load protein structures
2. Run binding site predictions
3. Access and analyze prediction results
4. Visualize binding sites with advanced plots
5. Compare different scoring methods
6. Analyze pocket properties

Let's get started!

## Setup

First, let's import the necessary modules and set up logging.
"""

import os
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import Image, display
from Bio.PDB import PDBParser, Selection

# Import ConSBind modules
from ConSBind.core.structure import ProteinStructure
from ConSBind.core.finder import ConsensusPocketFinder
from ConSBind.core.scoring import final_scoring
from ConSBind.output.output import save_predictions, save_pymol, save_chimera

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('ConSBind-Tutorial')

"""
## Load Protein Structure

First, we'll load a protein structure from a PDB file. ConSBind requires a valid PDB file as input.
For this tutorial, we'll use a sample PDB file from the data directory.

Let's define the path to our example PDB file:
"""

# Import ConSBind's file handling module
from ConSBind.input.file_handler import convert_ent_to_pdb, find_pdb_files
from pathlib import Path

# Create data directory if it doesn't exist
os.makedirs("data", exist_ok=True)

# Define the path to the PDB file
pdb_file = "data/pdb4dfr.pdb"
pdb_id = "4dfr"  # Dihydrofolate reductase

"""
## Downloading PDB Files

ConSBind can work with any valid PDB file. If you don't have a PDB file handy,
you can download one from the Protein Data Bank using BioPython's PDBList module.

Here's how to download a PDB file:
"""

# Check if the file exists
if not os.path.exists(pdb_file):
    print(f"PDB file {pdb_file} not found. Downloading it...")
    
    # Method 1: Using BioPython to download from PDB
    import Bio.PDB.PDBList
    pdblist = Bio.PDB.PDBList()
    
    # You can specify different file formats:
    # - "pdb" for PDB format (default)
    # - "mmCif" for mmCIF format
    # - "xml" for PDBML/XML format
    # - "mmtf" for MMTF format
    print(f"Downloading PDB ID: {pdb_id} in PDB format...")
    pdblist.retrieve_pdb_file(pdb_id, file_format="pdb", pdir="data")
    
    # PDBList downloads files with .ent extension
    ent_file = Path(f"data/pdb{pdb_id}.ent")
    
    if ent_file.exists():
        print(f"Downloaded PDB file to {ent_file}")
        # Convert .ent to .pdb using ConSBind's file handler
        pdb_file = str(convert_ent_to_pdb(ent_file))
        print(f"Converted to PDB format: {pdb_file}")
    else:
        # Try direct .pdb file
        pdb_file = f"data/{pdb_id}.pdb"
        if not os.path.exists(pdb_file):
            print(f"Could not find {ent_file} or {pdb_file}. Please download the PDB file manually.")
            sys.exit(1)
else:
    print(f"Using existing PDB file: {pdb_file}")

"""
## Loading the Protein Structure

Now that we have a PDB file, we can load it using ConSBind's ProteinStructure class.
This class handles parsing the PDB file and preparing it for binding site prediction.
"""

# Load the protein structure
try:
    print(f"Loading protein structure from {pdb_file}...")
    protein = ProteinStructure(pdb_file)
    print("Protein structure loaded successfully!")
    
    # Display basic information about the protein
    print(f"Protein ID: {protein.pdb_id}")
    print(f"Number of chains: {len(list(protein.structure.get_chains()))}")
    print(f"Number of residues: {len(list(protein.structure.get_residues()))}")
    print(f"Number of atoms: {len(list(protein.structure.get_atoms()))}")
except Exception as e:
    print(f"Error loading protein structure: {str(e)}")
    sys.exit(1)

"""
## Protein Structure Analysis

Before predicting binding sites, let's analyze the protein structure to understand its properties.
"""

def analyze_protein_structure(protein):
    """Analyze and visualize protein structure properties"""
    # Get basic protein information
    structure = protein.structure
    chains = list(structure.get_chains())
    residues = list(structure.get_residues())
    atoms = list(structure.get_atoms())
    
    print(f"Protein ID: {protein.pdb_id}")
    print(f"Number of chains: {len(chains)}")
    print(f"Number of residues: {len(residues)}")
    print(f"Number of atoms: {len(atoms)}")
    
    # Analyze residue composition
    residue_types = {}
    for residue in residues:
        res_name = residue.get_resname()
        if res_name in residue_types:
            residue_types[res_name] += 1
        else:
            residue_types[res_name] = 1
    
    # Exclude non-standard residues
    standard_residues = {
        'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine', 
        'ASP': 'Aspartic acid', 'CYS': 'Cysteine', 'GLN': 'Glutamine', 
        'GLU': 'Glutamic acid', 'GLY': 'Glycine', 'HIS': 'Histidine', 
        'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine', 
        'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline', 
        'SER': 'Serine', 'THR': 'Threonine', 'TRP': 'Tryptophan', 
        'TYR': 'Tyrosine', 'VAL': 'Valine'
    }
    
    # Filter for standard residues
    std_residue_counts = {k: v for k, v in residue_types.items() if k in standard_residues}
    
    # Plot residue composition
    plt.figure(figsize=(12, 6))
    plt.bar(std_residue_counts.keys(), std_residue_counts.values(), color='skyblue')
    plt.title('Amino Acid Composition')
    plt.xlabel('Amino Acid')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    
    # Categorize residues by property
    hydrophobic = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    polar = ['SER', 'THR', 'ASN', 'GLN', 'CYS']
    charged_positive = ['LYS', 'ARG', 'HIS']
    charged_negative = ['ASP', 'GLU']
    special = ['GLY', 'PRO']
    
    categories = {
        'Hydrophobic': sum(std_residue_counts.get(aa, 0) for aa in hydrophobic),
        'Polar': sum(std_residue_counts.get(aa, 0) for aa in polar),
        'Positive': sum(std_residue_counts.get(aa, 0) for aa in charged_positive),
        'Negative': sum(std_residue_counts.get(aa, 0) for aa in charged_negative),
        'Special': sum(std_residue_counts.get(aa, 0) for aa in special)
    }
    
    # Plot residue categories
    plt.figure(figsize=(10, 6))
    plt.pie(categories.values(), labels=categories.keys(), autopct='%1.1f%%', 
            shadow=True, startangle=140, colors=['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0'])
    plt.axis('equal')
    plt.title('Residue Properties Distribution')
    plt.tight_layout()
    plt.show()
    
    return residue_types

# Analyze the protein structure
residue_composition = analyze_protein_structure(protein)

"""
## Creating a Binding Site Prediction Function

Next, let's define a function to predict binding sites:
"""

def predict_binding_sites(protein, protein_type='enzyme', probe_radius=1.4, min_size=5):
    """
    Predict binding sites for a protein structure.
    
    Parameters:
    -----------
    protein : ProteinStructure
        The protein structure to analyze
    protein_type : str
        Type of protein ('enzyme', 'transporter', 'receptor', or 'unknown')
    probe_radius : float
        Probe radius for cavity detection
    min_size : int
        Minimum pocket size
        
    Returns:
    --------
    list
        List of pocket dictionaries with binding site information
    """
    print(f"Predicting binding sites for {protein.pdb_id}...")
    
    # Create a pocket finder
    finder = ConsensusPocketFinder()
    
    # Find pockets using geometric method
    geometric_pockets = finder.find_pockets_geometric(protein, probe_radius=probe_radius, min_size=min_size)
    print(f"Found {len(geometric_pockets)} pockets using geometric method")
    
    # Find pockets using energy-based method
    energy_pockets = finder.find_pockets_energy(protein)
    print(f"Found {len(energy_pockets)} pockets using energy-based method")
    
    # Combine pockets from different methods
    consensus_pockets = finder.combine_pockets(geometric_pockets, energy_pockets)
    print(f"Combined into {len(consensus_pockets)} consensus pockets")
    
    # Apply final scoring and filtering
    final_pockets = final_scoring(consensus_pockets, protein_type)
    print(f"Final selection: {len(final_pockets)} significant pockets")
    
    # Add residue information to each pocket
    for pocket in final_pockets:
        # Get residues within 8Å of pocket center
        pocket_residues = protein.get_pocket_residues(pocket, radius=8.0)
        
        # Format residues as chain:resname+resnum (e.g., "A:ALA123")
        formatted_residues = []
        for chain_id, res_id, res_name in pocket_residues:
            # Handle insertion codes if present
            if isinstance(res_id, tuple):
                res_num = res_id[1]  # Extract residue number from tuple
            else:
                res_num = res_id
            
            formatted_residues.append(f"{chain_id}:{res_name}{res_num}")
        
        pocket['residues'] = formatted_residues
    
    return final_pockets

# Run the prediction
pockets = predict_binding_sites(protein, protein_type='enzyme')

"""
## Visualizing Binding Site Predictions

Lets visualize the binding site predictions in various ways:
"""

def display_binding_sites(pockets, protein):
    """
    Display detailed information about predicted binding sites.
    """
    if not pockets:
        print("No binding sites found.")
        return {}
    
    print("\n" + "="*60)
    print(" " * 18 + "PREDICTED BINDING SITES")
    print("="*60)
    
    site_data = {
        'site_id': [],
        'consensus_score': [],
        'final_score': [],
        'knowledge_score': [],
        'size': [],
        'center_x': [],
        'center_y': [],
        'center_z': [],
        'methods': [],
        'residue_count': []
    }
    
    for i, pocket in enumerate(pockets):
        print(f"\nBinding Site {i+1} Details:")
        print("-" * 60 + "\n")
        
        # Extract scores
        binding_score = pocket.get('final_score', 0)
        consensus_score = pocket.get('consensus_score', 0)
        knowledge_score = pocket.get('knowledge_score', 0)
        
        print(f"Binding Potential Score: {binding_score:.2f}")
        print(f"Consensus Score: {consensus_score:.2f}")
        print(f"Knowledge-based Score: {knowledge_score:.2f}")
        
        # Extract size and center
        size = pocket.get('size', 0)
        center = pocket.get('center', [0, 0, 0])
        
        print(f"Size: {size}")
        print(f"Center: [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}]")
        
        # Extract detection methods
        methods = pocket.get('methods', [])
        print(f"Detection Methods: {', '.join(methods)}")
        
        # Extract residues
        residues = pocket.get('residues', [])
        print("\nBinding Site Residues:")
        for residue in residues:
            print(f"  {residue}")
        
        # Store data for visualization
        site_data['site_id'].append(i+1)
        site_data['consensus_score'].append(consensus_score)
        site_data['final_score'].append(binding_score)
        site_data['knowledge_score'].append(knowledge_score)
        site_data['size'].append(size)
        site_data['center_x'].append(center[0])
        site_data['center_y'].append(center[1])
        site_data['center_z'].append(center[2])
        site_data['methods'].append(', '.join(methods))
        site_data['residue_count'].append(len(residues))
    
    return site_data

# Display the binding site details
site_data = display_binding_sites(pockets, protein)

"""
## Advanced Visualizations

Let's create some advanced visualizations to better understand the binding site predictions:
"""

def create_advanced_visualizations(site_data, protein):
    """Create advanced visualizations of binding site data"""
    if not site_data or not site_data['site_id']:
        print("No binding site data available for visualization.")
        return
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(site_data)
    
    # 1. Score Comparison Bar Chart
    plt.figure(figsize=(12, 6))
    
    # Set width of bars
    barWidth = 0.25
    
    # Set positions of bars on X axis
    r1 = np.arange(len(df['site_id']))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    
    # Create bars
    plt.bar(r1, df['consensus_score'], width=barWidth, label='Consensus Score', color='skyblue')
    plt.bar(r2, df['final_score'], width=barWidth, label='Final Score', color='lightgreen')
    plt.bar(r3, df['knowledge_score'], width=barWidth, label='Knowledge Score', color='salmon')
    
    # Add labels and title
    plt.xlabel('Binding Site ID', fontweight='bold')
    plt.ylabel('Score Value', fontweight='bold')
    plt.title('Comparison of Different Scoring Methods')
    plt.xticks([r + barWidth for r in range(len(df['site_id']))], df['site_id'])
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # 2. 3D Scatter Plot of Binding Site Centers
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot with size proportional to binding site size
    scatter = ax.scatter(df['center_x'], df['center_y'], df['center_z'], 
                         c=df['final_score'], s=df['size']*2, 
                         cmap='viridis', alpha=0.7)
    
    # Add labels
    for i, txt in enumerate(df['site_id']):
        ax.text(df['center_x'][i], df['center_y'][i], df['center_z'][i], 
                f"Site {txt}", size=8)
    
    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Final Score')
    
    # Set labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Z coordinate')
    ax.set_title('3D Distribution of Binding Sites')
    
    # Set equal aspect ratio
    max_range = np.array([
        df['center_x'].max() - df['center_x'].min(),
        df['center_y'].max() - df['center_y'].min(),
        df['center_z'].max() - df['center_z'].min()
    ]).max() / 2.0
    
    mid_x = (df['center_x'].max() + df['center_x'].min()) * 0.5
    mid_y = (df['center_y'].max() + df['center_y'].min()) * 0.5
    mid_z = (df['center_z'].max() + df['center_z'].min()) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    plt.tight_layout()
    plt.show()
    
    # 3. Correlation Heatmap
    plt.figure(figsize=(10, 8))
    numeric_cols = ['consensus_score', 'final_score', 'knowledge_score', 'size', 'residue_count']
    corr_matrix = df[numeric_cols].corr()
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title('Correlation Between Binding Site Properties')
    plt.tight_layout()
    plt.show()
    
    # 4. Size vs. Score Scatter Plot with Regression Line
    plt.figure(figsize=(10, 6))
    sns.regplot(x='size', y='final_score', data=df, scatter_kws={'s': 80}, line_kws={'color': 'red'})
    plt.title('Relationship Between Binding Site Size and Final Score')
    plt.xlabel('Binding Site Size')
    plt.ylabel('Final Score')
    
    # Add site labels
    for i, txt in enumerate(df['site_id']):
        plt.annotate(f"Site {txt}", (df['size'][i], df['final_score'][i]), 
                    xytext=(5, 5), textcoords='offset points')
    
    plt.tight_layout()
    plt.show()
    
    # 5. Detection Methods Distribution
    # Extract method information
    method_counts = {'geometric': 0, 'energy': 0, 'both': 0}
    
    for method_str in df['methods']:
        methods = [m.strip() for m in method_str.split(',')]
        if 'geometric' in methods and 'energy' in methods:
            method_counts['both'] += 1
        elif 'geometric' in methods:
            method_counts['geometric'] += 1
        elif 'energy' in methods:
            method_counts['energy'] += 1
    
    # Create pie chart
    plt.figure(figsize=(8, 8))
    labels = ['Geometric Only', 'Energy Only', 'Both Methods']
    sizes = [method_counts['geometric'], method_counts['energy'], method_counts['both']]
    colors = ['#ff9999', '#66b3ff', '#99ff99']
    explode = (0.1, 0.1, 0.1)
    
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')
    plt.title('Distribution of Detection Methods')
    plt.tight_layout()
    plt.show()

# Create advanced visualizations
create_advanced_visualizations(site_data, protein)

"""
## Saving Prediction Results

Now, let's save the prediction results to files for further analysis and visualization:
"""

def save_results(pockets, protein, output_dir="results", generate_pymol=True, generate_chimera=True):
    """
    Save prediction results to files.
    
    Parameters:
    -----------
    pockets : list
        List of pocket dictionaries
    protein : ProteinStructure
        The protein structure
    output_dir : str
        Directory to save results
    generate_pymol : bool
        Whether to generate PyMOL visualization scripts
    generate_chimera : bool
        Whether to generate UCSF Chimera visualization scripts
        
    Returns:
    --------
    tuple
        Paths to the output files
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output prefix
    output_prefix = os.path.join(output_dir, protein.pdb_id)
    
    # Save predictions to text and PDB files
    output_file, output_pdb = save_predictions(pockets, protein, output_prefix)
    print(f"Saved prediction details to {output_file}")
    print(f"Saved modified PDB to {output_pdb}")
    
    # Generate visualization scripts based on user preferences
    pymol_script = None
    chimera_script = None
    
    if generate_pymol:
        pymol_script = save_pymol(pockets, protein, output_prefix, output_pdb)
        print(f"Saved PyMOL script to {pymol_script}")
    
    if generate_chimera:
        chimera_script = save_chimera(pockets, protein, output_prefix, output_pdb)
        print(f"Saved Chimera script to {chimera_script}")
    
    return output_file, output_pdb, pymol_script, chimera_script

"""
## Saving Results with Visualization Options

ConSBind can generate visualization scripts for both PyMOL and UCSF Chimera.
You can control which scripts are generated using the parameters.
"""

# Save the results with both PyMOL and Chimera scripts
output_files = save_results(pockets, protein, generate_pymol=True, generate_chimera=True)

print("\nVisualization instructions:")
print("1. To view in PyMOL: pymol results/4dfr/4dfr_pymol.pml")
print("2. To view in Chimera: chimera results/4dfr/4dfr_predicted.pdb results/4dfr/4dfr_chimera.cmd")

"""
## Analyzing Residue Properties in Binding Sites

Let's analyze the properties of residues in the binding sites:
"""

def analyze_binding_site_residues(pockets, protein):
    """
    Analyze the properties of residues in binding sites.
    
    Parameters:
    -----------
    pockets : list
        List of pocket dictionaries
    protein : ProteinStructure
        The protein structure
        
    Returns:
    --------
    dict
        Dictionary with residue analysis data
    """
    if not pockets:
        print("No binding sites to analyze.")
        return {}
    
    # Define residue property groups
    residue_properties = {
        'Hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
        'Polar': ['SER', 'THR', 'ASN', 'GLN', 'CYS'],
        'Positive': ['LYS', 'ARG', 'HIS'],
        'Negative': ['ASP', 'GLU'],
        'Special': ['GLY', 'PRO']
    }
    
    # Initialize data structure
    analysis_data = {
        'site_id': [],
        'total_residues': [],
        'hydrophobic': [],
        'polar': [],
        'positive': [],
        'negative': [],
        'special': []
    }
    
    # Analyze each pocket
    for i, pocket in enumerate(pockets):
        residues = pocket.get('residues', [])
        
        # Count residues by property
        property_counts = {prop: 0 for prop in residue_properties}
        
        for residue_str in residues:
            # Extract residue name (format: "A:ALA123")
            parts = residue_str.split(':')
            if len(parts) == 2:
                res_name = parts[1][:3]  # First 3 characters are the residue name
                
                # Count by property
                for prop, res_list in residue_properties.items():
                    if res_name in res_list:
                        property_counts[prop] += 1
                        break
        
        # Store data
        analysis_data['site_id'].append(i+1)
        analysis_data['total_residues'].append(len(residues))
        analysis_data['hydrophobic'].append(property_counts['Hydrophobic'])
        analysis_data['polar'].append(property_counts['Polar'])
        analysis_data['positive'].append(property_counts['Positive'])
        analysis_data['negative'].append(property_counts['Negative'])
        analysis_data['special'].append(property_counts['Special'])
    
    # Create DataFrame
    df = pd.DataFrame(analysis_data)
    
    # Calculate percentages
    for prop in ['hydrophobic', 'polar', 'positive', 'negative', 'special']:
        df[f'{prop}_pct'] = df[prop] / df['total_residues'] * 100
    
    # Plot residue composition for each binding site
    plt.figure(figsize=(14, 8))
    
    # Set width of bars
    barWidth = 0.15
    
    # Set positions of bars on X axis
    r1 = np.arange(len(df['site_id']))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]
    r5 = [x + barWidth for x in r4]
    
    # Create bars
    plt.bar(r1, df['hydrophobic_pct'], width=barWidth, label='Hydrophobic', color='#ff9999')
    plt.bar(r2, df['polar_pct'], width=barWidth, label='Polar', color='#66b3ff')
    plt.bar(r3, df['positive_pct'], width=barWidth, label='Positive', color='#99ff99')
    plt.bar(r4, df['negative_pct'], width=barWidth, label='Negative', color='#ffcc99')
    plt.bar(r5, df['special_pct'], width=barWidth, label='Special', color='#c2c2f0')
    
    # Add labels and title
    plt.xlabel('Binding Site ID', fontweight='bold')
    plt.ylabel('Percentage (%)', fontweight='bold')
    plt.title('Residue Property Distribution in Binding Sites')
    plt.xticks([r + barWidth*2 for r in range(len(df['site_id']))], df['site_id'])
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Create a heatmap of residue properties
    plt.figure(figsize=(12, 8))
    heatmap_data = df[['hydrophobic_pct', 'polar_pct', 'positive_pct', 'negative_pct', 'special_pct']]
    heatmap_data.index = [f"Site {i}" for i in df['site_id']]
    sns.heatmap(heatmap_data, annot=True, cmap='YlGnBu', fmt='.1f')
    plt.title('Residue Property Distribution Heatmap')
    plt.tight_layout()
    plt.show()
    
    return df

# Analyze binding site residues
residue_analysis = analyze_binding_site_residues(pockets, protein)

"""
## Comparing Multiple Proteins (Example)

In a real-world scenario, you might want to compare binding sites across multiple proteins.
Here's how you could set up such an analysis:
"""

def compare_multiple_proteins():
    """
    Demonstrates how to batch process multiple proteins with ConSBind.
    This function downloads example PDB files, processes each protein to predict binding sites,
    and compares the results across different protein types.
    """
    # Create a directory for multiple proteins
    multi_pdb_dir = "data/multiple_proteins"
    os.makedirs(multi_pdb_dir, exist_ok=True)
    
    # Define proteins to analyze (ID, description, type)
    proteins = [
        ("1HSG", "HIV-1 protease", "enzyme"),
        ("4DFR", "Dihydrofolate reductase", "enzyme"),
        ("3PTB", "Bovine trypsin", "enzyme"),
        ("1STP", "Streptavidin", "receptor")
    ]
    
    # Download proteins if needed
    for pdb_id, description, protein_type in proteins:
        pdb_file = os.path.join(multi_pdb_dir, f"{pdb_id}.pdb")
        if not os.path.exists(pdb_file):
            print(f"Downloading {pdb_id} - {description}")
            try:
                import Bio.PDB.PDBList
                pdblist = Bio.PDB.PDBList()
                pdblist.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=multi_pdb_dir)
                
                # Convert .ent to .pdb if needed
                ent_file = Path(f"{multi_pdb_dir}/pdb{pdb_id.lower()}.ent")
                if ent_file.exists():
                    from ConSBind.input.file_handler import convert_ent_to_pdb
                    pdb_file = convert_ent_to_pdb(ent_file)
            except Exception as e:
                print(f"Error downloading {pdb_id}: {str(e)}")
    
    # Find all PDB files in the directory
    from ConSBind.input.file_handler import find_pdb_files
    pdb_files = find_pdb_files(Path(multi_pdb_dir))
    
    if not pdb_files:
        print("No PDB files found in the directory.")
        return
    
    # Process each protein and collect results
    results = []
    
    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem
        print(f"Processing {pdb_id}...")
        
        # Find protein info
        protein_info = next((p for p in proteins if p[0] == pdb_id), (pdb_id, "Unknown", "unknown"))
        description = protein_info[1]
        protein_type = protein_info[2]
        
        try:
            # Load protein structure
            protein = ProteinStructure(str(pdb_file))
            
            # Predict binding sites
            pockets = predict_binding_sites(protein, protein_type=protein_type)
            
            if not pockets:
                print(f"No binding sites found for {pdb_id}")
                continue
                
            # Calculate metrics
            num_pockets = len(pockets)
            max_score = max([p.get('final_score', 0) for p in pockets])
            avg_size = sum([p.get('size', 0) for p in pockets]) / num_pockets if num_pockets > 0 else 0
            
            # Analyze residue properties
            residue_counts = {'hydrophobic': 0, 'polar': 0, 'charged': 0, 'other': 0}
            total_residues = 0
            
            # Define residue property groups
            residue_types = {
                'hydrophobic': ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
                'polar': ['SER', 'THR', 'ASN', 'GLN', 'CYS'],
                'charged': ['LYS', 'ARG', 'HIS', 'ASP', 'GLU']
            }
            
            # Print first pocket's residues for debugging
            if pockets and len(pockets) > 0 and 'residues' in pockets[0]:
                first_pocket = pockets[0]
                print(f"First pocket has {len(first_pocket.get('residues', []))} residues")
                if len(first_pocket.get('residues', [])) > 0:
                    print(f"Example residue: {first_pocket['residues'][0]}")
            
            for pocket in pockets:
                pocket_residues = pocket.get('residues', [])
                if not pocket_residues:
                    print(f"Warning: No residues found for pocket in {pdb_id}")
                    continue
                    
                total_residues += len(pocket_residues)
                
                for residue in pocket_residues:
                    # Extract residue name from the residue identifier
                    # Format is typically chain:resname+resnum (e.g., "A:ALA123")
                    try:
                        chain_id, res_info = residue.split(':')
                        res_name = res_info[:3]  # First 3 characters are the residue name
                        
                        if res_name in residue_types['hydrophobic']:
                            residue_counts['hydrophobic'] += 1
                        elif res_name in residue_types['polar']:
                            residue_counts['polar'] += 1
                        elif res_name in residue_types['charged']:
                            residue_counts['charged'] += 1
                        else:
                            residue_counts['other'] += 1
                    except Exception as e:
                        print(f"Error parsing residue {residue}: {str(e)}")
                        residue_counts['other'] += 1
            
            # Calculate percentages
            residue_percentages = {}
            if total_residues > 0:
                for res_type, count in residue_counts.items():
                    residue_percentages[f"{res_type}_percent"] = (count / total_residues) * 100
                print(f"Residue composition: {residue_counts} (total: {total_residues})")
            else:
                print("Warning: No residues found for analysis")
            
            # Store results
            results.append({
                'protein_id': pdb_id,
                'description': description,
                'protein_type': protein_type,
                'num_pockets': num_pockets,
                'max_score': max_score,
                'avg_pocket_size': avg_size,
                'total_residues': total_residues,
                **residue_percentages
            })
            
        except Exception as e:
            print(f"Error processing {pdb_id}: {str(e)}")
    
    # Create DataFrame with results
    if not results:
        print("No results collected.")
        return
        
    results_df = pd.DataFrame(results)
    
    # Create visualizations
    if len(results_df) > 1:
        # Binding sites by protein
        plt.figure(figsize=(10, 5))
        sns.barplot(x='protein_id', y='num_pockets', hue='protein_type', data=results_df)
        plt.title('Number of Binding Sites by Protein')
        plt.xlabel('Protein ID')
        plt.ylabel('Number of Binding Sites')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()
        
        # Pocket size vs. score
        plt.figure(figsize=(10, 5))
        sns.scatterplot(x='avg_pocket_size', y='max_score', hue='protein_type', 
                        size='num_pockets', sizes=(50, 200), data=results_df)
        
        # Add protein labels
        for i, row in results_df.iterrows():
            plt.annotate(row['protein_id'], 
                        (row['avg_pocket_size'], row['max_score']),
                        xytext=(5, 5), textcoords='offset points')
            
        plt.title('Relationship Between Pocket Size and Score')
        plt.xlabel('Average Pocket Size')
        plt.ylabel('Maximum Score')
        plt.tight_layout()
        plt.show()
        
        # Residue composition
        if 'hydrophobic_percent' in results_df.columns:
            plt.figure(figsize=(10, 5))
            # Reshape data for grouped bar chart
            property_data = pd.melt(results_df, 
                                  id_vars=['protein_id'],
                                  value_vars=['hydrophobic_percent', 'polar_percent', 'charged_percent'],
                                  var_name='property', value_name='percentage')
            
            # Create grouped bar chart
            sns.barplot(x='protein_id', y='percentage', hue='property', data=property_data)
            plt.title('Binding Site Residue Properties by Protein')
            plt.xlabel('Protein ID')
            plt.ylabel('Percentage of Residues')
            plt.xticks(rotation=45)
            plt.legend(title='Residue Property')
            plt.tight_layout()
            plt.show()
    
    return results_df

# Compare multiple proteins
comparison_results = compare_multiple_proteins()

"""
## Command-Line Usage

ConSBind can also be used from the command line, which is useful for batch processing multiple proteins or integrating into automated workflows. Let's explore how to use the command-line interface:

```bash
# Basic usage with a single PDB file
python main.py path/to/your/protein.pdb

# Specify protein type (improves prediction accuracy)
python main.py path/to/your/protein.pdb --protein_type enzyme

# Process all PDB files in a directory
python main.py path/to/pdb_directory/

# Customize output directory
python main.py path/to/your/protein.pdb --output_dir my_results

# Adjust prediction parameters
python main.py path/to/your/protein.pdb --min_size 3 --probe_radius 1.6 --grid_spacing 0.8

# Control visualization script generation
python main.py path/to/your/protein.pdb --generate_pymol --generate_chimera
```

### Command-Line Arguments

Here's a summary of the available command-line arguments:

| Argument | Description | Default |
|----------|-------------|--------|
| `input_path` | Input PDB file or directory containing PDB files | (Required) |
| `--output_dir` | Output directory | results |
| `--generate_pymol` | Generate PyMOL visualization script | True |
| `--generate_chimera` | Generate UCSF Chimera visualization script | True |
| `--min_size` | Minimum pocket size | 5 |
| `--probe_radius` | Probe radius for cavity detection | 1.4 |
| `--grid_spacing` | Grid spacing for energy calculations | 1.0 |
| `--consensus_threshold` | Minimum consensus score for reliable pockets | 1.5 |
| `--protein_type` | Type of protein (enzyme, transporter, receptor, unknown) | unknown |

### Example: Batch Processing Multiple Proteins

Here's a practical example of how to use ConSBind to batch process multiple proteins:

```python
import os
import subprocess
from pathlib import Path

# Directory containing PDB files
pdb_dir = "data/multiple_proteins"
os.makedirs(pdb_dir, exist_ok=True)

# List of interesting proteins to analyze
proteins = [
    ("1HSG", "HIV-1 protease (enzyme)"),
    ("4DFR", "Dihydrofolate reductase (enzyme)"),
    ("3PTB", "Bovine trypsin (enzyme)"),
    ("1STP", "Streptavidin (receptor)")
]

# Download proteins if needed
for pdb_id, description in proteins:
    pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    if not os.path.exists(pdb_file):
        print(f"Downloading {pdb_id} - {description}")
        import Bio.PDB.PDBList
        pdblist = Bio.PDB.PDBList()
        pdblist.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=pdb_dir)
        
        # Convert .ent to .pdb if needed
        ent_file = Path(f"{pdb_dir}/pdb{pdb_id.lower()}.ent")
        if ent_file.exists():
            from ConSBind.input.file_handler import convert_ent_to_pdb
            pdb_file = convert_ent_to_pdb(ent_file)
            print(f"Converted to PDB format: {pdb_file}")

# Run ConSBind on all proteins
output_dir = "results/multiple_comparison"
cmd = f"python main.py {pdb_dir} --output_dir {output_dir} --generate_pymol --generate_chimera"
print(f"Running: {cmd}")
# In a real script, you would use:
# subprocess.run(cmd, shell=True, check=True)
print("This would process all PDB files in the directory with both PyMOL and Chimera visualization scripts")
```

## Visualizing Results with PyMOL and Chimera

ConSBind generates visualization scripts for both PyMOL and UCSF Chimera. Here's how to use them:

### PyMOL Visualization

The PyMOL script colors the predicted binding sites and displays them as surfaces:

```bash
# Open PyMOL and load the script
pymol results/4dfr/4dfr_pymol.pml
```

The script automatically:
1. Loads the protein structure
2. Creates a surface representation
3. Colors each binding site with a different color
4. Labels the binding sites
5. Sets up a clean visualization with transparent surfaces

### UCSF Chimera Visualization

The Chimera script provides an alternative visualization:

```bash
# Open Chimera with the PDB file and run the script
chimera results/4dfr/4dfr_predicted.pdb results/4dfr/4dfr_chimera.cmd
```

The script:
1. Loads the protein structure
2. Creates a surface representation
3. Colors binding sites with distinct colors
4. Labels residues in each binding site
5. Sets up camera angles for optimal viewing

### Interpreting Visualization Results

When examining the visualizations, look for:

1. **Pocket location**: Are the binding sites in biologically relevant locations?
2. **Pocket size**: Larger pockets often indicate more significant binding sites
3. **Residue composition**: Check if known catalytic or binding residues are included
4. **Surface properties**: Examine the electrostatic and hydrophobic properties of the pocket

## Troubleshooting Common Issues

### Warning: Statistical Filtering Failed

As mentioned earlier, you might see this warning:
```
WARNING - Statistical filtering failed: list index out of range. Using all pockets.
```

This occurs when there aren't enough pockets to perform statistical filtering. Solutions:
- Try different protein types (--protein_type parameter)
- Adjust the minimum pocket size (--min_size parameter)
- Change the probe radius (--probe_radius parameter)

### Missing or Incorrect Binding Sites

If ConSBind doesn't find the expected binding sites:

1. **Check protein preparation**: Ensure the PDB file is complete and has proper hydrogen atoms
2. **Try different parameters**: Adjust probe_radius, min_size, and grid_spacing
3. **Consider protein type**: Different protein types use different scoring algorithms
4. **Examine all pockets**: Sometimes lower-ranked pockets may be biologically relevant

### Visualization Issues

If you encounter problems with visualization scripts:

1. **PyMOL version compatibility**: The scripts are designed for PyMOL 2.0+
2. **Chimera version compatibility**: The scripts are designed for UCSF Chimera 1.14+
3. **File paths**: Ensure the paths in the scripts match your actual file locations
4. **Manual loading**: You can always load the predicted PDB file manually and color by B-factor

## Conclusion

In this enhanced tutorial, we've demonstrated how to:

1. Load protein structures using ConSBind
2. Predict binding sites using consensus approaches
3. Analyze and visualize binding site properties
4. Save prediction results for external visualization
5. Compare properties across different binding sites
6. Set up a framework for multi-protein analysis
7. Use ConSBind from the command line
8. Process multiple proteins in batch mode
9. Visualize results with PyMOL and UCSF Chimera
10. Troubleshoot common issues

ConSBind provides a powerful framework for binding site prediction that combines geometric and energy-based approaches. By using the consensus approach, you can achieve more reliable predictions than with any single method.

These techniques can be applied to your own protein structures to identify potential binding sites
for drug discovery, protein engineering, or functional analysis.

For more information, refer to the ConSBind documentation and the original research papers.
"""

print("Tutorial completed successfully!")
