@echo off
conda env list | findstr /c:"ngsolve" >nul 2>&1
if %errorlevel% neq 0 (
    echo conda creating new enviroment ...
    conda create -n ngsolve python=3.8
    call conda activate ngsolve
) else (
    REM Package is already installed
    echo enviroment is already created
    call conda activate ngsolve
)

REM Define the list of packages to install
set "packages=numpy scipy matplotlib pytz datetime pyevtk ngsolve"

REM Loop through each package
for %%P in (%packages%) do (
    REM Check if the package is already installed
    python -c "import %%P" >nul 2>&1 && echo %%P is already installed || (
        echo Installing %%P...
        pip install %%P
    )
)

echo All packages have been processed.

set "PYTHONPATH=%PYTHONPATH%;.\Packages"
