#!/bin/bash
# filepath: deploy.sh

# MLegS Package Deployment Script
# Automatically sets up MLegS on supercomputers or local machines

set -e  # Exit on any error

echo "========================================="
echo "MLegS Package Deployment Script"
echo "========================================="

# Function to detect system type
detect_system() {
    if command -v srun >/dev/null 2>&1 || command -v sbatch >/dev/null 2>&1 || command -v qsub >/dev/null 2>&1; then
        echo "supercomputer"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo "macos"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "linux"
    else
        echo "unknown"
    fi
}

# Function to load modules on supercomputer
load_supercomputer_modules() {
    echo "Detected supercomputer environment"
    echo "Loading required modules..."
    
    # Check if module command exists
    if command -v module >/dev/null 2>&1; then
        module load intel impi conda 2>/dev/null || {
            echo "Warning: Could not load some modules. Trying alternatives..."
            module load intel/2021.4 2>/dev/null || echo "Intel module not found"
            module load impi/2021.4 2>/dev/null || echo "Intel MPI module not found"
            module load anaconda3 2>/dev/null || module load miniconda3 2>/dev/null || echo "Conda module not found"
        }
        echo "Loaded modules:"
        module list 2>&1 | head -10
    else
        echo "Warning: Module system not available"
    fi
}

# Function to build lib_fm
build_lib_fm() {
    echo "========================================="
    echo "Building lib_fm library..."
    echo "========================================="
    
    if [ ! -d "lib_fm" ]; then
        echo "Error: lib_fm directory not found!"
        exit 1
    fi
    
    cd lib_fm
    echo "Cleaning previous build..."
    make clean || echo "Warning: make clean failed"
    
    echo "Building lib_fm..."
    if make all; then
        echo "lib_fm built successfully"
    else
        echo "Error: lib_fm build failed!"
        exit 1
    fi
    
    cd ..
}

# Function to initialize main package
init_package() {
    echo "========================================="
    echo "Initializing MLegS package..."
    echo "========================================="
    
    if [ ! -f "Makefile" ]; then
        echo "Error: Makefile not found in root directory!"
        exit 1
    fi
    
    echo "Running make init..."
    if make init; then
        echo "Package initialized successfully"
    else
        echo "Error: Package initialization failed!"
        exit 1
    fi
}

# Function to setup conda environment
# Function to setup conda environment using environment.yml
setup_conda_env() {
    echo "========================================="
    echo "Setting up Python environment..."
    echo "========================================="
    
    # Check if conda is available
    if ! command -v conda >/dev/null 2>&1; then
        echo "Error: Conda not found! Please install Anaconda or Miniconda first."
        echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
        return 1
    fi
    
    # Check for environment.yml
    if [ ! -f "environment.yml" ]; then
        echo "Error: environment.yml not found in repository!"
        echo "Please create an environment.yml file in your repository root."
        exit 1
    fi
    
    echo "Creating conda environment from environment.yml..."
    
    # Remove existing environment if it exists
    conda env remove -n mlegs -y 2>/dev/null || true
    
    # Create environment from yml file
    if conda env create -f environment.yml; then
        echo "Conda environment 'mlegs' created successfully from environment.yml"
    else
        echo "Error: Failed to create conda environment from environment.yml"
        return 1
    fi
    
    echo "Activating environment and setting up Jupyter extensions..."
    
    # Activate environment 
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate mlegs
    
    # Install additional Jupyter extensions
    echo "Installing Jupyter extensions..."
    jupyter nbextension enable --py widgetsnbextension --sys-prefix 2>/dev/null || true
    
    echo "Python environment setup complete!"
    echo ""
    echo "To activate the environment, run:"
    echo "  conda activate mlegs"
    echo ""
    echo "To start Jupyter Lab, run:"
    echo "  jupyter lab"
    echo ""
    echo "To update environment from environment.yml:"
    echo "  conda env update -f environment.yml"
    
    conda deactivate
}

# Main deployment function
main() {
    echo "Starting MLegS deployment..."
    
    # Check if we're in the right directory
    if [ ! -f "run_ivp.sh" ] && [ ! -d "src" ]; then
        echo "Error: This doesn't appear to be the MLegS root directory!"
        echo "Please run this script from the MLegS package root."
        exit 1
    fi
    
    # Detect system type
    SYSTEM_TYPE=$(detect_system)
    echo "Detected system type: $SYSTEM_TYPE"
    
    # Load modules on supercomputer
    if [ "$SYSTEM_TYPE" = "supercomputer" ]; then
        load_supercomputer_modules
    fi
    
    # Build lib_fm
    build_lib_fm
    
    # Initialize package
    init_package
    
    # Setup Python environment
    echo ""
    read -p "Do you want to set up the Python analysis environment? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        setup_conda_env
    fi
    
    echo ""
    echo "========================================="
    echo "Deployment completed successfully!"
    echo "========================================="
    echo ""
}

# Run main function
main "$@"