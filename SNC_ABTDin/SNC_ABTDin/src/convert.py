import os
import pandas as pd


#with os.scandir('C:\Users\marcr\source\repos\SNC_ABTDin\data\graph') as entries:
   # for entry in entries:
df = pd.read_csv(r'C:\Users\marcr\source\repos\SNC_ABTDin\data\graph\Slashdot0902.txt', sep=',', header=None)
df.to_csv(r'C:\Users\marcr\source\repos\SNC_ABTDin\data\graph\modified\1Slashdot0902.txt', sep=' ', header=None, index=None)
