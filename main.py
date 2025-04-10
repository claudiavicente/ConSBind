#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ConSBind - Consensus Structural Binding site predictor
=====================================================
This program predicts potential binding sites in protein structures
using a consensus of geometric and energy-based approaches.

Input: 
    - Single PDB file (--pdb option)
    - Directory containing multiple PDB files (--pdb_dir option)

Output: 
    - Text file listing amino acids involved in each detected site
    - Modified PDB file for visualization in molecular graphics software
    - PyMOL script for visualization (optional, enabled by default)
    - UCSF Chimera script for visualization (optional, enabled by default)

Output Directory Structure:
    - For single PDB file:
      results/pdb_basename/pdb_basename_predictions.txt
      results/pdb_basename/pdb_basename_predicted.pdb
      results/pdb_basename/pdb_basename_pymol.pml
      results/pdb_basename/pdb_basename_chimera.py (if --generate_chimera is used)
      
    - For directory of PDB files:
      results/dir_basename/pdb1_basename/pdb1_basename_predictions.txt
      results/dir_basename/pdb1_basename/pdb1_basename_predicted.pdb
      results/dir_basename/pdb1_basename/pdb1_basename_pymol.pml
      results/dir_basename/pdb1_basename/pdb1_basename_chimera.py (if --generate_chimera is used)
      results/dir_basename/pdb2_basename/...
 """

import os
import sys
import logging
import argparse
from pathlib import Path
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from colorama import Fore, Style, init

from ConSBind.core.structure import ProteinStructure
from ConSBind.core.finder import ConsensusPocketFinder
from ConSBind.core.scoring import final_scoring
from ConSBind.output.output import save_predictions, save_pymol, save_chimera
from ConSBind.input.file_handler import detect_input_type, find_pdb_files, create_output_path

# Initialize colorama for cross-platform colored terminal output
init(autoreset=True)

# Configure logging
class ColoredFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': Fore.BLUE,
        'INFO': Fore.GREEN,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRITICAL': Fore.RED + Style.BRIGHT
    }
    
    def format(self, record):
        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{Style.RESET_ALL}"
            if levelname == 'WARNING':
                record.msg = f"{Fore.YELLOW}{record.msg}{Style.RESET_ALL}"
            elif levelname in ['ERROR', 'CRITICAL']:
                record.msg = f"{Fore.RED}{record.msg}{Style.RESET_ALL}"
        return super().format(record)

# Custom handler for tqdm compatibility
class TqdmLoggingHandler(logging.StreamHandler):
    def __init__(self):
        super().__init__()
    
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

# Configure tqdm to work with logging
tqdm.set_lock(tqdm.get_lock())

# Set up the logger
logger = logging.getLogger('ConSBind')
logger.setLevel(logging.INFO)
logger.propagate = False  # Prevent propagation to root logger

# Remove existing handlers
for handler in logger.handlers:
    logger.removeHandler(handler)

# Create custom tqdm-compatible handler with colored formatter
tqdm_handler = TqdmLoggingHandler()
tqdm_handler.setLevel(logging.INFO)
formatter = ColoredFormatter(
    fmt=f'{Fore.LIGHTBLACK_EX}%(asctime)s{Style.RESET_ALL} - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
tqdm_handler.setFormatter(formatter)
logger.addHandler(tqdm_handler)

# Configure other loggers to prevent duplicate messages
logging.getLogger('tqdm').setLevel(logging.WARNING)  # Reduce tqdm log noise

def process_single_pdb(pdb_file, output_path, args):
    """
    Process a single PDB file for binding site prediction
    
    Parameters:
    -----------
    pdb_file : str
        Path to the PDB file
    output_path : str
        Path for output files
    args : argparse.Namespace
        Command line arguments
    
    Returns:
    --------
    bool
        Success or failure
    """
    try:
        pdb_basename = os.path.basename(pdb_file)
        
        # Create output directory if it doesn't exist
        os.makedirs(output_path, exist_ok=True)
        
        # Set output prefix
        output_prefix = os.path.join(output_path, os.path.splitext(pdb_basename)[0])
        
        # Print initial information before starting progress bar
        logger.info(f"Processing: {Fore.CYAN}{pdb_basename}{Style.RESET_ALL}")
        
        # Create a progress bar for the processing steps
        steps = [
            "Loading protein structure",
            "Finding geometric pockets",
            "Finding energy-based pockets",
            "Combining results",
            "Scoring pockets",
            "Saving predictions"
        ]
        
        if args.generate_pymol:
            steps.append("Generating PyMOL visualization")
        
        if args.generate_chimera:
            steps.append("Generating Chimera visualization")
        
        # Initialize progress bar with position=0 to keep it at the bottom
        with tqdm(total=len(steps), desc=f"Processing {pdb_basename}", 
                  bar_format="{l_bar}{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                  position=0, leave=True, dynamic_ncols=True, 
                  file=sys.stdout) as pbar:
            
                # Load protein structure
                pbar.set_description(f"Loading {pdb_basename}")
                protein = ProteinStructure(pdb_file)
                pbar.update(1)
                
                # Create consensus pocket finder
                pocket_finder = ConsensusPocketFinder()
                
                # Find pockets using geometric approach
                pbar.set_description(f"Finding geometric pockets in {pdb_basename}")
                geometric_pockets = pocket_finder.find_pockets_geometric(
                    protein, 
                    probe_radius=args.probe_radius,
                    min_size=args.min_size
                )
                pbar.update(1)
                
                # Find pockets using energy-based approach
                pbar.set_description(f"Finding energy-based pockets in {pdb_basename}")
                energy_pockets = pocket_finder.find_pockets_energy(
                    protein,
                    grid_spacing=args.grid_spacing
                )
                pbar.update(1)
                
                # Combine results
                pbar.set_description(f"Combining results for {pdb_basename}")
                consensus_pockets = pocket_finder.combine_pockets(protein, geometric_pockets, energy_pockets)
                pbar.update(1)
                
                # Adjust scores based on protein function
                pbar.set_description(f"Scoring pockets for {pdb_basename}")
                consensus_pockets = final_scoring(consensus_pockets, args.protein_type)
                
                # Ensure every pocket has a final score for proper sorting 
                for pocket in consensus_pockets:
                    if 'final_score' not in pocket:
                        pocket['final_score'] = pocket['consensus_score'] * 3.0 + pocket.get('score', 0) * 1.0
                        
                    # Ensure methods list exists
                    if 'methods' not in pocket:
                        pocket['methods'] = [pocket.get('method', 'unknown')]
                pbar.update(1)
                
                if not consensus_pockets:
                    logger.warning(f"No binding sites found for {pdb_basename}")
                else:
                    # Save predictions
                    pbar.set_description(f"Saving predictions for {pdb_basename}")
                    output_file, output_pdb = save_predictions(consensus_pockets, protein, output_prefix)
                    pbar.update(1)
                    
                    # Generate visualization scripts
                    if args.generate_pymol:
                        pbar.set_description(f"Generating PyMOL script for {pdb_basename}")
                        pymol_script = save_pymol(consensus_pockets, protein, output_prefix, output_pdb)
                        pbar.update(1)
                        
                    if args.generate_chimera:
                        pbar.set_description(f"Generating Chimera script for {pdb_basename}")
                        chimera_script = save_chimera(consensus_pockets, protein, output_prefix, output_pdb)
                        pbar.update(1)
        
        # Summary of results
        if consensus_pockets:
            num_pockets = len(consensus_pockets)
            logger.info(f"Found {Fore.YELLOW}{num_pockets}{Style.RESET_ALL} binding sites in {Fore.CYAN}{pdb_basename}{Style.RESET_ALL}")
        logger.info(f"Results saved to: {Fore.BLUE}{output_path}{Style.RESET_ALL}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error processing {pdb_basename}: {str(e)}")
        return False

def main():
    """Main function for binding site prediction"""
    parser = argparse.ArgumentParser(
        description='Predict protein binding sites using consensus approach',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    """
    )
    
    # Input argument - automatically detects file or directory
    parser.add_argument('input_path', help='Input PDB file or directory containing PDB files')
    
    # Output options
    parser.add_argument('--output_dir', default='results', 
                        help='Output directory (default: results)')
    parser.add_argument('--generate_pymol', action='store_true', default=False,
                        help='Generate PyMOL visualization script (default: False)')
    parser.add_argument('--generate_chimera', action='store_true', default=False,
                        help='Generate UCSF Chimera visualization script (default: False)')
    
    # Prediction parameters
    predict_group = parser.add_argument_group('Prediction Parameters')
    predict_group.add_argument('--min_size', type=int, default=5, 
                        help='Minimum pocket size (default: 5)')
    predict_group.add_argument('--probe_radius', type=float, default=1.4, 
                        help='Probe radius for cavity detection (default: 1.4)')
    predict_group.add_argument('--grid_spacing', type=float, default=1.0, 
                        help='Grid spacing for energy calculations (default: 1.0)')
    predict_group.add_argument('--consensus_threshold', type=float, default=1.5, 
                        help='Minimum consensus score for reliable pockets (default: 1.5)')
    predict_group.add_argument('--protein_type', 
                        choices=['enzyme', 'transporter', 'receptor', 'unknown'],
                        default='unknown', 
                        help='Type of protein for specialized detection (default: unknown)')
    
    args = parser.parse_args()
    
    # Create the base results directory
    base_output_dir = Path(args.output_dir)
    
    try:
        # Automatically detect if input is a file or directory
        input_type, input_path = detect_input_type(args.input_path)
        
        if input_type == 'file':
            # Process a single PDB file
            file_basename = input_path.name
            logger.info(f"Input: {Fore.CYAN}{file_basename}{Style.RESET_ALL}")
            
            # Create output directory structure
            output_path = create_output_path(input_path, base_output_dir)
            
            success = process_single_pdb(str(input_path), output_path, args)
            if not success:
                sys.exit(1)
                
        elif input_type == 'directory':
            # Process all PDB files in a directory
            dir_basename = input_path.name
            logger.info(f"Input directory: {Fore.CYAN}{dir_basename}{Style.RESET_ALL}")
            
            # Create output directory structure
            output_base_path = create_output_path(input_path, base_output_dir)
            
            # Find all PDB files in the directory
            pdb_files = find_pdb_files(input_path)
            
            if not pdb_files:
                logger.error(f"No PDB files found in directory: {dir_basename}")
                sys.exit(1)
            
            # Display summary of files to be processed
            logger.info(f"Found {Fore.YELLOW}{len(pdb_files)}{Style.RESET_ALL} PDB files to process")
            
            # Process each PDB file with a master progress bar
            success_count = 0
            # Process each PDB file with a master progress bar
            success_count = 0
            with tqdm(total=len(pdb_files), desc=f"Overall progress", 
                     bar_format="{l_bar}{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                     position=0, leave=True, dynamic_ncols=True,
                     file=sys.stdout) as master_pbar:
                
                    for pdb_file in pdb_files:
                        # Create output directory structure for each PDB file
                        pdb_output_path = output_base_path / pdb_file.stem
                        
                        # Process the file
                        if process_single_pdb(str(pdb_file), pdb_output_path, args):
                            success_count += 1
                        
                        # Update the master progress bar
                        master_pbar.update(1)
            
            # Final summary
            logger.info(f"{Fore.GREEN}Successfully processed {Fore.YELLOW}{success_count}{Fore.GREEN} out of {Fore.YELLOW}{len(pdb_files)}{Fore.GREEN} PDB files{Style.RESET_ALL}")
            logger.info(f"Results saved to: {Fore.BLUE}{output_base_path}{Style.RESET_ALL}")
            
            # Exit with error if no files were processed successfully
            if success_count == 0:
                logger.error("Failed to process any PDB files")
                sys.exit(1)
    
    except (FileNotFoundError, ValueError) as e:
        logger.error(str(e))
        sys.exit(1)

if __name__ == "__main__":
    main() 