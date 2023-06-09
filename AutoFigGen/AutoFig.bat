set "PYTHONPATH=%PYTHONPATH%;..\Packages"

@REM Example1: Standard-soliton in dimensional 1, convergence test

python Std1d_Spat_Conv.py
python Std1d_Temp_Conv.py
python Std1d_CPU_dt.py
python Std1d_ME_Conserv.py
python Std1d_En-Mass_Cp.py
python Std1d_Sol.py

@REM Example2: Bi-soliton long time behavior

python BiSoli1d_ShortTimeSol.py
python BiSoli1d_LongTimeSol.py
python BiSoli1d_LongT_ME_Err_Compare.py

@REM Example3: soliton in dimension 2

python 2dSoli_Spat_Conv.py
python 2dSoli_Temp_Conv.py
python 2dSoli_ME_Conserv.py