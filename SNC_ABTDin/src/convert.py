import os
import pandas as pd


#with os.scandir('C:\Users\marcr\source\repos\SNC_ABTDin\data\graph') as entries:
   # for entry in entries:
df = pd.read_csv(r'C:\Users\marcr\source\repos\SNC_ABTDin\data\graph\tvshow_edges.txt', sep=',', header=None)
df.to_csv(r'C:\Users\marcr\source\repos\SNC_ABTDin\data\graph\modified\tvshow_edges2.txt', sep=' ', header=None, index=None)
