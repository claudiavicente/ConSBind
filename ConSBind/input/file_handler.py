#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
File handling module for ConSBind
This module provides functions to detect and process input files and directories
"""

import os
import glob
import shutil
import tempfile
import logging
from pathlib import Path

# Use the same logger as main to prevent duplicate messages
logger = logging.getLogger('ConSBind')

def convert_ent_to_pdb(ent_file):
    """
    Convert .ent file to .pdb format
    
    Parameters:
    -----------
    ent_file : Path
        Path to the .ent file
    
    Returns:
    --------
    Path
        Path to the converted .pdb file (in a temporary directory if needed)
    """
    # If it's already a .pdb file, return it as is
    if ent_file.suffix.lower() == '.pdb':
        return ent_file
    
    logger.info(f"Converting .ent file to .pdb format: {os.path.basename(ent_file)}")
    
    # Create a temporary directory if needed
    temp_dir = Path(tempfile.gettempdir()) / "consbind_temp"
    temp_dir.mkdir(exist_ok=True)
    
    # Create a new .pdb file with the same content
    pdb_file = temp_dir / f"{ent_file.stem}.pdb"
    
    # Copy the content (no conversion needed, just rename)
    shutil.copy2(ent_file, pdb_file)
    
    logger.info(f"Converted to: {os.path.basename(pdb_file)}")
    return pdb_file

def detect_input_type(input_path):
    """
    Automatically detect if the input is a file or directory
    
    Parameters:
    -----------
    input_path : str
        Path to the input file or directory
    
    Returns:
    --------
    tuple
        (input_type, path)
        input_type: 'file' or 'directory'
        path: Path object of the input
    """
    path = Path(input_path)
    
    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")
    
    if path.is_file():
        # Check if it's a PDB file
        if path.suffix.lower() in ['.pdb', '.ent']:
            # Convert .ent to .pdb if needed
            if path.suffix.lower() == '.ent':
                path = convert_ent_to_pdb(path)
            return 'file', path
        else:
            raise ValueError(f"Input file is not a PDB file: {input_path}")
    
    elif path.is_dir():
        return 'directory', path
    
    else:
        raise ValueError(f"Input path is neither a file nor a directory: {input_path}")

def find_pdb_files(directory):
    """
    Find all PDB files in a directory and convert .ent files to .pdb format
    
    Parameters:
    -----------
    directory : str or Path
        Directory to search for PDB files
    
    Returns:
    --------
    list
        List of Path objects for PDB files (converted to .pdb if needed)
    """
    directory = Path(directory)
    
    if not directory.exists() or not directory.is_dir():
        raise ValueError(f"Invalid directory: {directory}")
    
    pdb_files = []
    
    # Find all files with PDB extensions
    for ext in ['.pdb', '.ent', '.PDB', '.ENT']:
        pdb_files.extend(directory.glob(f'*{ext}'))
    
    if not pdb_files:
        logger.warning(f"No PDB files found in directory: {directory}")
    else:
        logger.info(f"Found {len(pdb_files)} PDB files in {directory}")
    
    # Convert any .ent files to .pdb format
    converted_files = []
    for file in pdb_files:
        if file.suffix.lower() == '.ent':
            converted_files.append(convert_ent_to_pdb(file))
        else:
            converted_files.append(file)
    
    return converted_files

def create_output_path(input_path, base_output_dir):
    """
    Create appropriate output path based on input path
    
    Parameters:
    -----------
    input_path : Path
        Path to input file or directory
    base_output_dir : Path
        Base output directory
    
    Returns:
    --------
    Path
        Output path for results
    """
    if input_path.is_file():
        # For a single file: results/pdb_basename/
        pdb_basename = input_path.stem
        return base_output_dir / pdb_basename
    
    elif input_path.is_dir():
        # For a directory: results/dir_basename/
        dir_basename = input_path.name
        return base_output_dir / dir_basename
    
    else:
        raise ValueError(f"Invalid input path: {input_path}")
