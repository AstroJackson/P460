import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, os

os.chdir("/Users/jackson/Desktop/P460/Lab 1")

data = pd.read_excel("P460Lab1.xlsx", sheet_name="Data")
D, sigma_D = data["D"], data["dD"]
s01, sigma_s01 = data["S01"], data["dS01"]
s02, sigma_s02 = data["S02"], data["dS02"]
b1, sigma_b1 = data["b1"], data["db1"]
b2, sigma_b2 = data["b2"], data["db2"]

plt.plot(D, s01)
plt.show()
