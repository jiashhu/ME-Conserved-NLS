# 8192 - 87
# 4096 - 44            17:10
# 2048 - 28            10:59
# 1024 - 12:25         5:31
# 512  - 6:40          2:44
# 256  - 3:24          1:38
# 128  - 1:36          51
# 64   - 32            31
# 32   - 28            17
# 16   - 11            10
# 8    - 10            6

# tau = 1/(2**(np.array([3,4,5,6,7,8,9,10,11,12])))
# Direct_CPU_set = [6,10,17,31,51,98,164,
#                331,659,1030]

# 'Proj-Method'
#        1800 cpu      5000 cpu
# 8192 - 89
# 4096 - 41:55         17:14
# 2048 - 32:30         11:21
# 1024 - 16:14         5:34
# 512  - 8:10          49+120+4
# 256  - 4:30          1:43
# 128  - 2:32          54
# 64   - 1:30          33
# 32   - 42            18
# 16   - 30            11
# 8    - 8             7

# tau = 1/(2**(np.array([3,4,5,6,7,8,9,10,11,12])))
# Our_CPU_set = [7,11,18,33,54,103,173,
#                334,681,1034]

import pandas as pd
# pd.to_datetime().astype(int)/ 10**9
TimeProposedMethod = [
       ['2023-05-13 18:51:23', '2023-05-13 18:51:29', 8],
       ['2023-05-13 18:51:30', '2023-05-13 18:51:49', 32],
       ['2023-05-13 18:51:49', '2023-05-13 18:52:40', 128],
       ['2023-05-13 18:52:40', '2023-05-13 18:55:27', 512],
       ['2023-05-13 18:55:27', '2023-05-13 19:06:38', 2048],
       ['2023-05-14 16:12:49', '2023-05-14 16:13:01', 16],
       ['2023-05-14 16:13:02', '2023-05-14 16:13:32', 64],
       ['2023-05-14 16:13:33', '2023-05-14 16:15:17', 256],
       ['2023-05-14 16:15:18', '2023-05-14 16:20:53', 1024],
       ['2023-05-14 16:20:54', '2023-05-14 16:38:09', 4096]
]
df = pd.DataFrame(TimeProposedMethod, columns=['Btime','Etime','Tsteps'])
df['Btime'] = pd.to_datetime(df['Btime']).astype(int)/10**9
df['Etime'] = pd.to_datetime(df['Etime']).astype(int)/10**9
df['Ptime'] = df['Etime']-df['Btime']
df['tau'] = 1/df['Tsteps']

PP_CPU = [0.1130828857421875, 
              0.4572722911834717,
              1.4711434841156006,
              6.351077079772949,
              24.52998375892639,
              0.26062631607055664, 
              0.722142219543457,
              3.074068069458008,
              12.293888092041016,
              49.70548748970032,]
 
df['PP_CPU'] = PP_CPU

df.to_csv('/home/jiashhu/Collocation_NLS/ReviseCode/cputime.csv',index=False)
def TimeOfProjection():
       '''
              (NgsolveNew) jiashhu@AMAs3 Experiment-18-11-2022-BGN-Dziuk $ python3 /home/jiashhu/Collocation_NLS/ReviseCode/NLS_Collo_1d_Time_Conv.py
              Parameter of this example:
              dim = 1, n_collo = 2, order = 3, N_thres = 1e-10, ref_order = 1
              T = 1, N = 1024
              2023-05-13 18:51:23 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-13 18:51:24 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-13 18:51:25 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-13 18:51:25 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-13 18:51:26 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-13 18:51:27 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-13 18:51:28 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-13 18:51:29 +0800 Asia/Shanghai-Finished 70.0 per cent
              PPtime is 0.1130828857421875
              (0.031883895386665126+0j)
              8
              2023-05-13 18:51:30 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-13 18:51:31 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-13 18:51:33 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-13 18:51:35 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-13 18:51:37 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-13 18:51:39 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-13 18:51:41 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-13 18:51:43 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-13 18:51:45 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-13 18:51:47 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-13 18:51:49 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 0.4572722911834717
              (0.000291972911862182+0j)
              32
              2023-05-13 18:51:49 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-13 18:51:54 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-13 18:52:00 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-13 18:52:05 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-13 18:52:10 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-13 18:52:14 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-13 18:52:19 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-13 18:52:24 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-13 18:52:30 +0800 Asia/Shanghai-Finished 80.0 per cent
              no more timer available, reusing last one
              2023-05-13 18:52:35 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-13 18:52:40 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 1.4711434841156006
              (6.878981892510233e-06+0j)
              128
              2023-05-13 18:52:40 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-13 18:52:56 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-13 18:53:13 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-13 18:53:30 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-13 18:53:46 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-13 18:54:03 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-13 18:54:21 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-13 18:54:37 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-13 18:54:53 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-13 18:55:10 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-13 18:55:27 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 6.351077079772949
              (6.871323940047293e-06+0j)
              512
              2023-05-13 18:55:27 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-13 18:56:36 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-13 18:57:44 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-13 18:58:52 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-13 18:59:58 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-13 19:01:05 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-13 19:02:11 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-13 19:03:17 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-13 19:04:25 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-13 19:05:32 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-13 19:06:38 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 24.52998375892639
              (6.87290332473679e-06+0j)
              2048
              (NgsolveNew) jiashhu@AMAs3 Experiment-18-11-2022-BGN-Dziuk $ python3 /home/jiashhu/Collocation_NLS/ReviseCode/NLS_Collo_1d_Time_Conv.py
              Parameter of this example:
              dim = 1, n_collo = 2, order = 3, N_thres = 1e-10, ref_order = 1
              T = 1, N = 1024
              2023-05-14 16:12:49 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-14 16:12:50 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-14 16:12:51 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-14 16:12:52 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-14 16:12:54 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-14 16:12:54 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-14 16:12:56 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-14 16:12:57 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-14 16:12:58 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-14 16:13:00 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-14 16:13:01 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 0.26062631607055664
              (0.0029191451675724287+0j)
              16
              2023-05-14 16:13:02 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-14 16:13:05 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-14 16:13:07 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-14 16:13:11 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-14 16:13:14 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-14 16:13:17 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-14 16:13:20 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-14 16:13:23 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-14 16:13:26 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-14 16:13:29 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-14 16:13:32 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 0.722142219543457
              (3.209631017068662e-05+0j)
              64
              2023-05-14 16:13:33 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-14 16:13:43 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-14 16:13:53 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-14 16:14:03 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-14 16:14:14 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-14 16:14:25 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-14 16:14:35 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-14 16:14:46 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-14 16:14:56 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-14 16:15:06 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-14 16:15:17 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 3.074068069458008
              (6.862944447833157e-06+0j)
              256
              2023-05-14 16:15:18 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-14 16:15:52 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-14 16:16:25 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-14 16:16:59 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-14 16:17:32 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-14 16:18:06 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-14 16:18:40 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-14 16:19:13 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-14 16:19:46 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-14 16:20:20 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-14 16:20:53 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 12.293888092041016
              (6.872417050744516e-06+0j)
              1024
              2023-05-14 16:20:54 +0800 Asia/Shanghai-Finished 0.0 per cent
              2023-05-14 16:22:36 +0800 Asia/Shanghai-Finished 10.0 per cent
              2023-05-14 16:24:21 +0800 Asia/Shanghai-Finished 20.0 per cent
              2023-05-14 16:26:04 +0800 Asia/Shanghai-Finished 30.0 per cent
              2023-05-14 16:27:49 +0800 Asia/Shanghai-Finished 40.0 per cent
              2023-05-14 16:29:33 +0800 Asia/Shanghai-Finished 50.0 per cent
              2023-05-14 16:31:16 +0800 Asia/Shanghai-Finished 60.0 per cent
              2023-05-14 16:32:58 +0800 Asia/Shanghai-Finished 70.0 per cent
              2023-05-14 16:34:43 +0800 Asia/Shanghai-Finished 80.0 per cent
              2023-05-14 16:36:26 +0800 Asia/Shanghai-Finished 90.0 per cent
              2023-05-14 16:38:09 +0800 Asia/Shanghai-Finished 100.0 per cent
              PPtime is 49.70548748970032
              (6.873443896313806e-06+0j)
              4096
              NumSol saved Successful!
       '''
       pass
