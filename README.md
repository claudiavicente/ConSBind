# ConSBind: Consensus Structural Binding Site Predictor

ConSBind is a Python package for predicting protein binding sites using a consensus of geometric and energy-based approaches. The tool identifies potential binding pockets in protein structures and scores them based on multiple complementary criteria.

## Features

- Detection of protein binding sites using multiple methods
- Consensus scoring to improve prediction reliability
- Knowledge-based filters to identify catalytic sites and other functional regions
- Ability to customize prediction parameters for different protein types
- Output in standard PDB format with visualization support for PyMOL

## Installation

```bash
# Clone the repository
git clone https://github.com/sbiteam/ConSBind.git
cd ConSBind

# Install the package
pip install -e .
```

## Dependencies

- Python 3.6+
- NumPy
- SciPy
- BioPython
- scikit-learn
- DSSP (mkdssp binary needs to be in PATH)

### Installing DSSP

DSSP is required for surface calculations. Install it using one of the following methods:

```bash
# On Linux
sudo apt-get install dssp

# On macOS
brew install brewsci/bio/dssp

# Using conda (any platform)
conda install -c salilab dssp
```

Alternatively, you can build DSSP from source:
https://github.com/PDB-REDO/dssp

After installation, ensure the `mkdssp` binary is in your PATH or create a symlink called `mkdssp` that points to the DSSP executable:

```bash
# Example if the DSSP binary is named 'dssp'
sudo ln -s $(which dssp) /usr/local/bin/mkdssp
```

## Usage

### Basic Usage

```bash
consbind --pdb protein.pdb
```

This will generate:
- `protein_predictions.txt` - List of predicted binding sites and residues
- `protein_predicted.pdb` - Modified PDB file with binding sites
- `protein_pymol.pml` - PyMOL script for visualization

### Advanced Options

```bash
consbind --pdb protein.pdb --protein_type enzyme --max_pockets 3 --score_threshold 2.5
```

Available options:
- `--output` - Output file prefix (default: PDB filename)
- `--min_size` - Minimum pocket size (default: 5)
- `--probe_radius` - Probe radius for cavity detection (default: 1.4)
- `--grid_spacing` - Grid spacing for energy calculations (default: 1.0)
- `--max_pockets` - Maximum number of pockets to report (default: 5)
- `--score_threshold` - Minimum score threshold (default: 3.0)
- `--consensus_threshold` - Minimum consensus score (default: 1.5)
- `--protein_type` - Type of protein: enzyme, transporter, receptor, or unknown (default: unknown)

## Visualization

To visualize the results in PyMOL:

```bash
pymol protein_pymol.pml
```

## Binding Site Scores

ConSBind uses several score metrics:

1. **Consensus Score**: Measures agreement between different detection methods
2. **Binding Potential Score**: Composite score combining consensus and pocket characteristics
3. **Druggability Score**: Estimates the pocket's suitability for small molecule binding
4. **Knowledge-based Score**: Based on known binding site residue patterns

## Methodology

ConSBind combines multiple approaches to binding site detection:

1. **Geometric**: Identifies cavities and pockets in the protein structure
2. **Energy-based**: Evaluates the physicochemical properties of potential binding regions
3. **Knowledge-based**: Applies known patterns of binding site composition and architecture

The final predictions represent a consensus of these approaches, ranked by a combined score that prioritizes pockets detected by multiple methods.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
