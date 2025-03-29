#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Output and Visualization Module
==============================
This module handles the generation of output files and visualization.
"""

import os
import time
import numpy as np
import logging
from Bio.PDB import PDBIO

# Use the same logger as main to prevent duplicate messages
logger = logging.getLogger('ConSBind')

# Define colors based on score ranges for visualization
COLORS = {
    0.9: "forest",        # Excellent score (dark green)
    0.75: "lime",         # Very high score (light green)
    0.6: "limon",         # High score (yellow-green)
    0.45: "yellow",       # Medium score (yellow)
    0.3: "orange",        # Low-medium score (orange)
    0.15: "salmon",       # Low score (light red)
    0.0: "red"            # Very low score (red)
}


def save_predictions(pockets, protein, output_prefix):
    """Save predictions to text file and modified PDB with binding site indicators.

    Parameters
    ----------
    pockets : list
        List of pocket dictionaries containing binding site information
    protein : Protein object
        Protein object containing structure information
    output_prefix : str
        Prefix for output files

    Returns
    -------
    tuple
        Paths to the output text file and PDB file
    """
    output_file = f"{output_prefix}_predictions.txt"
    output_pdb = f"{output_prefix}_predicted.pdb"

    # Write text summary
    with open(output_file, 'w') as f:
        f.write("Predicted Binding Sites\n")
        f.write("======================\n\n")

        for i, pocket in enumerate(pockets, 1):
            # Get residues for this pocket
            residues = protein.get_pocket_residues(pocket)

            f.write(f"Site {i}:\n")

            # Show all methods that detected this site
            if 'methods' in pocket and len(pocket['methods']) > 0:
                f.write(f"Detection Methods: {', '.join(pocket['methods'])}\n")
            else:
                f.write(f"Detection Method: {pocket['method']}\n")

            f.write(f"Consensus Score: {pocket['consensus_score']:.2f}\n")
            f.write(f"Binding Potential Score: {pocket['final_score']:.2f}\n")

            if 'druggability' in pocket:
                f.write(f"Druggability: {pocket['druggability']:.2f}\n")
            if 'knowledge_score' in pocket:
                f.write(f"Knowledge-based Score: {pocket['knowledge_score']:.2f}\n")

            f.write(f"Size: {pocket['size']}\n")
            f.write(f"Center: {pocket['center'][0]:.3f}, {pocket['center'][1]:.3f}, {pocket['center'][2]:.3f}\n")

            f.write("\nBinding Site Residues:\n")
            for chain, resid, resname in protein.get_pocket_residues(pocket):
                f.write(f"  {chain}:{resname}{resid}\n")
            f.write("\n")

    # Get the structure data from the original PDB file
    pdb_io = PDBIO()
    pdb_io.set_structure(protein.structure)

    # Read the original PDB file
    with open(protein.pdb_file, 'r') as f:
        pdb_lines = f.readlines()

    # Create a modified PDB file with binding site indicators
    with open(output_pdb, 'w') as f:
        # Write all lines up to the END record (if it exists)
        end_found = False
        for line in pdb_lines:
            if line.startswith('END'):
                end_found = True
                break
            f.write(line)

        # Add REMARK records for binding sites
        f.write("\nREMARK   Predicted Binding Sites\n")
        f.write(f"REMARK   Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

        # Create a new chain identifier for binding sites
        binding_site_chain = 'X'  # Chain X for binding sites

        # Add dummy atoms for each prediction
        atom_num = 10000  # Start from a high number to avoid conflicts
        for i, pocket in enumerate(pockets, 1):
            f.write(f"REMARK   Site {i} - Method: {pocket['method']}, Consensus: {pocket['consensus_score']:.2f}, Binding Potential: {pocket['final_score']:.2f}\n")

            # Add the center point as a larger sphere - use ATOM instead of HETATM
            x, y, z = pocket['center']
            f.write(f"ATOM  {atom_num:5d}  O   SIT {binding_site_chain}{i:3d}    "
                   f"{x:8.3f}{y:8.3f}{z:8.3f}"
                   f"  1.00{pocket['consensus_score']:6.2f}           O\n")
            atom_num += 1

            # Add smaller spheres for sample points in the cluster if available
            points = pocket.get('points', [])
            if isinstance(points, np.ndarray) and len(points) > 0:
                # Select a subset of points for visualization (max 20)
                if len(points) > 20:
                    indices = np.random.choice(len(points), 20, replace=False)
                    points = points[indices]

                for j, point in enumerate(points):
                    x, y, z = point
                    f.write(f"ATOM  {atom_num:5d}  H   SIT {binding_site_chain}{i:3d}    "
                           f"{x:8.3f}{y:8.3f}{z:8.3f}"
                           f"  1.00  0.00           H\n")
                    atom_num += 1

        # Add connectivity and END records
        f.write("TER\n")
        f.write("END\n")

    logger.info(f"Predictions saved to {output_file}")
    logger.info(f"Modified PDB saved to {output_pdb}")

    return output_file, output_pdb


def save_pymol(pockets, protein, output_prefix, output_pdb=None):
    """Generate a PyMOL script for visualizing binding sites.

    Parameters
    ----------
    pockets : list
        List of pocket dictionaries containing binding site information
    protein : Protein object
        Protein object containing structure information
    output_prefix : str
        Prefix for output files
    output_pdb : str, optional
        Path to the predicted PDB file. If None, will be generated from output_prefix

    Returns
    -------
    str
        Path to the PyMOL script file
    """
    # If output_pdb is not provided, construct it from the prefix
    if output_pdb is None:
        output_pdb = f"{output_prefix}_predicted.pdb"

    # Generate PyMOL script for visualization
    pymol_script = f"{output_prefix}_pymol.pml"
    with open(pymol_script, 'w') as f:
        # Get the base name of the PDB file to use as the object name
        pdb_name = os.path.splitext(os.path.basename(output_pdb))[0]

        f.write(f"# PyMOL script for visualizing predicted binding sites\n")
        # Use just the filename instead of the full path for better portability
        pdb_filename = os.path.basename(output_pdb)
        f.write(f"load {pdb_filename}, main_obj\n")  # Explicitly name the object
        f.write(f"hide everything\n")

        # Set background color to white and adjust display settings
        f.write(f"bg_color white\n")
        f.write(f"set antialias, 2\n")  # Better antialiasing
        f.write(f"set light_count, 2\n")  # Add more lights for better shading
        f.write(f"set specular, 0.1\n")  # Reduce specular highlights
        f.write(f"set sphere_quality, 2\n")  # Higher quality spheres
        f.write(f"set cartoon_fancy_helices, 1\n")  # Better looking helices
        f.write(f"set depth_cue, 0\n")  # No depth cueing with white background

        # Label settings
        f.write(f"set label_size, 10\n")  # Global label size
        f.write(f"set label_position, (0,0,0)\n")  # Center positioning
        f.write(f"set label_color, black\n")  # Black text
        f.write(f"set label_depth_mask, 0\n")  # Force labels in front regardless of depth
        f.write(f"set float_labels, on\n")  # Make labels float in front of objects

        # Show protein cartoon but not binding site indicators
        f.write(f"show cartoon, main_obj and not chain X\n")
        f.write(f"color gray80, main_obj and not chain X\n")
        f.write(f"set cartoon_transparency, 0.5\n\n")

        # Show ligands if present (but not binding site indicators)
        f.write("# Show ligands in magenta\n")
        f.write(f"select ligands, main_obj and hetatm and not resn HOH and not chain X\n")
        f.write("show sticks, ligands\n")
        f.write("color magenta, ligands\n\n")

        # Add representation for each binding site
        for i, pocket in enumerate(pockets, 1):
            # Get score and determine color - normalize the consensus_score for color selection
            normalized_score = min(1.0, max(0.0, pocket['consensus_score'] / 5.0))
            site_color = "red"  # default color
            for threshold, color in COLORS.items():
                if normalized_score >= threshold:
                    site_color = color
                    break

            f.write(f"# Binding site {i}\n")
            f.write(f"select site_{i}_center, (main_obj and chain X and resi {i} and name O)\n")
            f.write(f"select site_{i}_points, (main_obj and chain X and resi {i} and name H)\n")

            # Show spheres for center (black) and points (colored by score)
            f.write(f"show spheres, site_{i}_center\n")
            f.write(f"color black, site_{i}_center\n")  # Center point always black
            f.write(f"set sphere_scale, 1.0, site_{i}_center\n")  # Make center points larger

            # Only display points if they exist
            f.write(f"show spheres, site_{i}_points\n")
            f.write(f"color {site_color}, site_{i}_points\n")  # Points colored by score
            f.write(f"set sphere_scale, 0.6, site_{i}_points\n")  # Make cluster points smaller

            # Select and display residues
            residues = protein.get_pocket_residues(pocket)
            if residues:
                residue_sel = " or ".join([f"(main_obj and chain {chain} and resi {resid})" for chain, resid, _ in residues])
                f.write(f"select site_{i}_res, ({residue_sel})\n")
                f.write(f"show sticks, site_{i}_res\n")
                f.write(f"color {site_color}, site_{i}_res\n")

                # Create surface for the binding site
                f.write(f"create site_{i}_surface, site_{i}_res\n")
                f.write(f"show surface, site_{i}_surface\n")
                f.write(f"set surface_color, {site_color}, site_{i}_surface\n")
                f.write(f"set transparency, 0.3, site_{i}_surface\n\n")

            # Add label at the center of the binding site
            x, y, z = pocket['center']

            # Create a small empty sphere specifically for labels
            f.write(f"# Add label sphere at binding site center\n")
            f.write(f"pseudoatom label_obj_{i}, pos=[{x}, {y}, {z}]\n")
            f.write(f"set sphere_scale, 0.01, label_obj_{i}\n")  # Make it tiny but still present
            f.write(f"color black, label_obj_{i}\n")
            f.write(f"label label_obj_{i}, \"Site {i} (Score: {pocket['final_score']:.1f})\"\n\n")

            # Add site details as comments
            f.write(f"# Site {i} details:\n")
            f.write(f"# Consensus Score: {pocket['consensus_score']:.2f}\n")
            f.write(f"# Binding Potential Score: {pocket['final_score']:.2f}\n")
            f.write(f"# Size: {pocket['size']}\n")

            # Show all detection methods in comments
            if 'methods' in pocket and len(pocket['methods']) > 0:
                f.write(f"# Detection Methods: {', '.join(pocket['methods'])}\n")
            else:
                f.write(f"# Method: {pocket['method']}\n")

            f.write("\n")

        # Final view settings
        f.write("# Set up view\n")
        f.write("orient\n")
        f.write(f"zoom main_obj, 20\n")
        f.write(f"center main_obj\n\n")

        # Ray trace settings
        f.write("# Ray trace settings for better quality\n")
        f.write("set ray_shadows, 1\n")  # Enable shadows with white background
        f.write("set ray_shadow_decay_factor, 0.1\n")  # Softer shadows
        f.write("set ray_trace_mode, 1\n")
        f.write("set ray_trace_color, black\n")  # Black outlines
        f.write("set ray_trace_gain, 0.8\n")  # Slightly lower gain for white background
        f.write("set ray_opaque_background, on\n\n")  # Ensure solid white background in images

        # Group all objects
        f.write("# Group all objects for easier manipulation\n")
        f.write(f"group binding_sites_centers, site_*_center\n")
        f.write(f"group binding_sites_points, site_*_points\n")
        f.write(f"group binding_sites_residues, site_*_res\n")
        f.write(f"group binding_sites, binding_sites_*\n")
        f.write(f"group labels, label_obj_*\n")
        f.write(f"group surfaces, *_surface\n\n")

        # Save session
        f.write("# Save session\n")
        f.write(f"save {pdb_name}_binding_sites.pse\n")

    logger.info(f"PyMOL script saved to {pymol_script}")
    return pymol_script


def save_chimera(pockets, protein, output_prefix, output_pdb=None):
    """Generate a UCSF Chimera script for visualizing binding sites.

    Parameters
    ----------
    pockets : list
        List of pocket dictionaries containing binding site information
    protein : Protein object
        Protein object containing structure information
    output_prefix : str
        Prefix for output files
    output_pdb : str, optional
        Path to the predicted PDB file. If None, will be generated from output_prefix

    Returns
    -------
    str
        Path to the Chimera script file
    """
    # If output_pdb is not provided, construct it from the prefix
    if output_pdb is None:
        output_pdb = f"{output_prefix}_predicted.pdb"

    # Create a BILD file for the spheres
    bild_file = f"{output_prefix}_spheres.bild"
    with open(bild_file, 'w') as bf:
        bf.write(".transparency 0.0\n")  # Make spheres solid

        # Add a sphere for each binding site
        for i, pocket in enumerate(pockets, 1):
            # Get score and determine color
            score = pocket['final_score']

            # Simple color scheme with RGB values (required for BILD format)
            if score > 4.0:
                r, g, b = 0.13, 0.55, 0.13  # forest green
            elif score > 3.0:
                r, g, b = 0.0, 1.0, 0.0  # green
            elif score > 2.0:
                r, g, b = 1.0, 1.0, 0.0  # yellow
            elif score > 1.0:
                r, g, b = 1.0, 0.65, 0.0  # orange
            elif score > 0.5:
                r, g, b = 1.0, 0.27, 0.0  # orange red
            else:
                r, g, b = 1.0, 0.0, 0.0  # red

            # Get center coordinates
            x, y, z = pocket['center']

            # Add comments for binding site information
            bf.write(f".comment Binding site {i} (Score: {score:.2f})\n")

            # Set color and draw sphere
            bf.write(f".color {r:.2f} {g:.2f} {b:.2f}\n")
            bf.write(f".sphere {x:.2f} {y:.2f} {z:.2f} 2.0\n\n")

    # Create a simple Chimera script file to open the BILD file
    chimera_script = f"{output_prefix}_chimera.cmd"
    with open(chimera_script, 'w') as f:
        # Write header and instructions
        f.write("# Chimera script for visualizing binding sites\n")
        f.write("# To use this script:\n")
        f.write("# 1. Open Chimera\n")
        f.write("# 2. Open your PDB file\n")
        f.write("# 3. Go to Tools -> General Controls -> Command Line\n")
        f.write("# 4. Type: open /path/to/this/file.cmd\n\n")

        # Basic display settings
        f.write("# Basic display settings\n")
        f.write("background solid white\n")

        # Open the BILD file with the spheres
        bild_basename = os.path.basename(bild_file)
        f.write(f"# Open the BILD file with binding site spheres\n")
        f.write(f"open {bild_basename}\n\n")

        # Final view settings
        f.write("# Center the view\n")
        f.write("focus\n")

    logger.info(f"UCSF Chimera BILD file saved to {bild_file}")
    logger.info(f"UCSF Chimera script saved to {chimera_script}")
    return chimera_script
