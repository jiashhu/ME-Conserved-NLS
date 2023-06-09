export PYTHONPATH=$PYTHONPATH:../Packages

## Running experiments for testing spatial convergence of FEM degree 1,2,3
## And generate .npy files collecting H1 errs 
# python NLS_Collo_1d.py 1d_Standard_Soliton/Spat_Conv/Spat_Conv_o1.json
# python NLS_Collo_1d.py 1d_Standard_Soliton/Spat_Conv/Spat_Conv_o2.json
# python NLS_Collo_1d.py 1d_Standard_Soliton/Spat_Conv/Spat_Conv_o3.json

## Read all .npy files according to the parameters .json files and collect the infH1 error into a Sum.. file
# python Data_Analyse_H1err.py 1d_Standard_Soliton/Spat_Conv

## Running experiments for testing temporal convergence of collocation degree 2,3,4
## And generate .npy files collecting H1 errs 
# python NLS_Collo_1d.py 1d_Standard_Soliton/Temp_Conv/Temp_Conv_o2.json
# python NLS_Collo_1d.py 1d_Standard_Soliton/Temp_Conv/Temp_Conv_o3.json
# python NLS_Collo_1d.py 1d_Standard_Soliton/Temp_Conv/Temp_Conv_o4.json

## Read all .npy files according to the parameters .json files and collect the infH1 error into a Sum.. file
# python Data_Analyse_H1err.py 1d_Standard_Soliton/Temp_Conv

# python NLS_Collo_1d.py 1d_Standard_Soliton/ME_Conserv/ME_Conserv.json
# python NLS_Collo_1d.py 1d_Standard_Soliton/En-Mass_Cp/En-Mass_PPC.json
# python NLS_Collo_1d.py 1d_Standard_Soliton/En-Mass_Cp/En-Mass_GL.json

# python NLS_Collo_2d.py 2d_Soliton/Spat_Conv/Spat_Conv_o1.json
# python NLS_Collo_2d.py 2d_Soliton/Spat_Conv/Spat_Conv_o2.json
# python NLS_Collo_2d.py 2d_Soliton/Spat_Conv/Spat_Conv_o3.json

# python Data_Analyse_H1err.py 2d_Soliton/Spat_Conv

python NLS_Collo_2d.py 2d_Soliton/Temp_Conv/Temp_Conv_o2.json
python NLS_Collo_2d.py 2d_Soliton/Temp_Conv/Temp_Conv_o3.json
python NLS_Collo_2d.py 2d_Soliton/Temp_Conv/Temp_Conv_o4.json

# python Data_Analyse_H1err.py 2d_Soliton/Temp_Conv

