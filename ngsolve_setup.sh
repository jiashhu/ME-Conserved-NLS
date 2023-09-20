#!/bin/bash

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Installing Miniconda..."
    
    # Download the Miniconda installation script
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    
    # Make the script executable
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    
    # Run the installation script
    ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
    
    # Add Conda to the system PATH
    export PATH="$HOME/miniconda/bin:$PATH"
    
    # Activate Conda
    source $HOME/miniconda/bin/activate
    
    echo "Miniconda installed and activated."
else
    echo "Conda found. Activating Conda environment..."
    
    # Activate Conda
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate
    
    echo "Conda environment activated."
fi

if conda env list | grep -q "ngsolve"; then
    # Environment exists
    echo "Environment already created."
    source activate ngsolve
else
    # Environment does not exist
    echo "Creating new environment..."
    conda create -n ngsolve python=3.8
    source activate ngsolve
fi

#!/bin/bash

# Define the list of packages to install
packages=("numpy" "scipy" "matplotlib" "pytz" "datetime" "pyevtk" "ngsolve")

# Loop through each package
for package in "${packages[@]}"; do
    # Check if the package is already installed
    if python -c "import $package" > /dev/null 2>&1; then
        echo "$package is already installed."
    else
        echo "Installing $package..."
        pip install "$package"
    fi
done

echo All packages have been processed.

export "PYTHONPATH=%PYTHONPATH%:./Packages"
