from __future__ import print_function
import mdtraj as md
import numpy as np
import itertools

t=md.load('1yrc_added.pdb')                         # reference

# Extract interface index
group1 = range(0,85)
group2 = range(85,98)
pairs = list(itertools.product(group1, group2))
A=md.compute_contacts(t, pairs, scheme='closest-heavy')

H=[] 
A1=min(A[0])   # distance array
A2=A[1]        # residue pairs 

for x in range(len(A1)):
	if A1[x]<=1.0: # condition 1
		H.append(x)

# Extract the index from A2 which satisfy condition1 and save it to A3
A3=A2[H]
A4=[]
A5=[]

# Extract the index of residues which belong to protein
for x in range(len(A3)):
	A4.append(A3[x][0])

A5=list(set(A4)) # the protein residues' index which consist of the interface

# Appending peptide residue index
for x in group2: 
	A5.append(x)

#Extract the reference frame based on crystal structure
atom_to_keep = t.topology.select('backbone')
t.atom_slice(atom_to_keep, inplace=True)
print(t.topology.atoms)
A6=[]
for item in A5:
	A6.append(str(item))
A7=['resid '+ x for x in A6]
A8=[]

for x in A7:
	m=t.topology.select(x)
	n=m.tolist()
	A8=A8+n

print(A8) # all the interface residues' atoms 
atom_to_keep = A8
#print(A8)
t.atom_slice(atom_to_keep, inplace=True)
t.save('reference_interface.h5') # the reference frame only alpha carbon

##############################################simulation result###################################################
t=md.load('trajectory.00.dcd', top='topol.prmtop')   # trajectory
#Extract the alpha carbon of complex

atom_to_keep = t.topology.select('backbone')
t.atom_slice(atom_to_keep,inplace=True) # this acts must be done on currently loaded trajectory

atom_to_keep = A8
t.atom_slice(atom_to_keep, inplace=True)
t.save('interface.h5') # the modified trajectory

#load the interface trajectory and calculate the rmsd
t1=md.load('interface.h5')
t2=md.load('reference_interface.h5')

t1.center_coordinates # trajectory
t2.center_coordinates # reference

t2.superpose(t1)
rmsd_interface = md.rmsd(t1,t2, 0,precentered=True) #using the first frame of t1 as reference

#plot the results
from matplotlib.pylab import *
figure()
plot(t1.time, rmsd_interface,'r', label='interface rmsd')
plt.show()
   


