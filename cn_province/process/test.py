import numpy as np

f = open("cn_province.txt","r")
string = f.read().split()
print int(string[3])
f.close()
