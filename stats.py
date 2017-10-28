import pandas
import pickle
import numpy as np

with open('data/randomResults/HLdata.txt', 'rb') as f:
	x = pickle.load(f)
with open('data/randomResults/SM_E.txt', 'rb') as f:
	y = pickle.load(f)

X = x.values()
Y = y.values()

print np.corrcoef(X,Y)