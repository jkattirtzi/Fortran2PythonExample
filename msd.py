import numpy
import msd4py
Ifile=open("ACNpcom_nocell.xyz",'r')
natoms=43 # the number of atoms will read in, not necessarily total number in a frame
n=0
i=0
Txarray=[]
Tyarray=[]
Tzarray=[]
for line in Ifile:
    l=line.split()
    if l[0]=="X":
        n+=1
        x=float(l[1])
        y=float(l[2])
        z=float(l[3])
        Txarray.append(x)
        Tyarray.append(y)
        Tzarray.append(z)
        if n==natoms:
            i+=1
            n=0
Txarray=numpy.array(Txarray)
Tyarray=numpy.array(Tyarray)
Tzarray=numpy.array(Tzarray)
msd4py.msd.calc_msd(Txarray,Tyarray,Tzarray,natoms,i)
