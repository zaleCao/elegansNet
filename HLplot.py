import pickle
import numpy as np
import matplotlib.pyplot as plt

with open('data/randomResults/HLdata.txt', 'rb') as f:
	HL = pickle.load(f)	

x = HL.keys()
y = HL.values()

fig, ax = plt.subplots()
line1, = ax.plot(x, y, linewidth=1,)
plt.axis([0, 20, -0.20, 0.20])

fig.suptitle('Hub - Non-Hub Value ', fontsize=14)
plt.ylabel('HN Value')
plt.xlabel('Time Delay h')

plt.show()