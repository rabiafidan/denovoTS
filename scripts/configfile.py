#! /usr/bin/python

import sys
import os
path=sys.argv[1]
os.chdir(path)

#%%
c=[]
with open("children.txt", "r") as child:
        for line in child:
            c.append(line.strip())
# %%
m=[]
with open("mother.txt", "r") as mot:
        for line in mot:
            m.append(line.strip())
# %%
f=[]
with open("father.txt", "r") as fat: 
        for line in fat:
            f.append(line.strip())

# %%
with open("../config.yaml", "w") as con:
    con.write("trios:\n")
    for i in range(1,len(f)+1):
        con.write("        '"+str(i)+"\':")
        con.write("\n")
        con.write("                mother: "+m[i-1])
        con.write("\n")
        con.write("                father: "+f[i-1])
        con.write("\n")
        con.write("                child: "+c[i-1])
        con.write("\n")
        con.write("\n")
# %%
