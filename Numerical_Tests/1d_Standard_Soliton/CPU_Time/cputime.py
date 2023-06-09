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
       49.70548748970032]
 
df['PP_CPU'] = PP_CPU

df.to_csv('cputime.csv',index=False)
