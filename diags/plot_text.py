import matplotlib.pyplot as plt
import sys


#Figure 1
plt.figure()
file = open(sys.argv[1], "r")
xs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ys = []
for line in file:
 ys.append(float(line.split()[-1]))
plt.plot(xs, ys)
title = sys.argv[1].split("/")[-1].replace(".txt","")
plt.title(title)
plt.savefig("transp1.ps")

#Figure 2
plt.figure()
file = open(sys.argv[2], "r")
xs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ys = []
for line in file:
 ys.append(float(line.split()[-1]))
plt.plot(xs, ys)
title = sys.argv[2].split("/")[-1].replace(".txt","")
plt.title(title)
plt.savefig("transp2.ps")

#Figure 3
plt.figure()
file = open(sys.argv[3], "r")
xs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
ys = []
for line in file:
 ys.append(float(line.split()[-1]))
plt.plot(xs, ys)
title = sys.argv[3].split("/")[-1].replace(".txt","")
plt.title(title)
plt.savefig("transp2.ps")

plt.show()
