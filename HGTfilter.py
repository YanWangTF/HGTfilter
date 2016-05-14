#! /usr/bin/python

import os
import sys
import numpy
import math
import getopt

def usage():
	print """
HGTfilter.py: reads a blast output and builds two dictionaries with best Host_Hit and 
best Fungi_Hit of each gene; E-values of the Host_Hit and Fungi_Hit from the same gene 
will be compared to sort out the potential HGT elements. Three files will be produced: 
absolute_HGT (the genes have only Host_Hit), FineFilted_HGT (the absolute value of 
logEvalue of the Host_Hit is at least two times bigger than the Fungi_Hit), and 
CoarseFilted_HGT (genes of the Host_Hit have better Evalue than Fungi_Hit's)

Usage: HGTfilter.py [-h] <Blastoutputname>

-h                  print this help message

<Blastoutputname>   the file has to be the format 6; the database sequence header has to 
					start with "Host_" or "Fungi_"

"""

o, a = getopt.getopt(sys.argv[1:],'-output:h')
opts = {}

for k,v in o:
	opts[k]=v
if '-h' in opts.keys():
	usage(); sys.exit()

gene_list = []
coarseHGT=[]
fineHGT=[]
abHGT=[]
query_tag=[0]
base_tag=[0]

Harp_Host={}
HarpHost_E={}
Harp_Fungi={}
HarpFungi_E={}

Blastoutput= sys.argv[1];

try:
	t=open(Blastoutput)
except IOError:
	print("File %s does not exit!!!" % Blastoutput)

filecontent=open(Blastoutput)

for line in filecontent:
	gene= line.split()[0]
	hit= line.split()[1]
	evalue= line.split()[10]

	if gene !=  query_tag[-1]:
		query_tag.append(gene)
		base_tag.append(hit.split('_')[0])
		if hit.split('_')[0] =='Host':
			coarseHGT.append(gene)
			Harp_Host[gene] = hit
			HarpHost_E[gene] = evalue
		else:
			Harp_Fungi[gene] = hit
			HarpFungi_E[gene] = evalue

	else:
		if gene in list(Harp_Host.keys()) and gene in list(Harp_Fungi.keys()):
			pass
		else:
			if hit.split('_')[0] ==base_tag[-1]:
				pass
			else:
				base_tag.append(hit.split('_')[0])
				if hit.split('_')[0] =='Host':
					Harp_Host[gene] = hit
					HarpHost_E[gene] = evalue
				else:
					Harp_Fungi[gene] = hit
					HarpFungi_E[gene] = evalue


for m in list(HarpHost_E.keys()):
	if HarpHost_E[m]=='0.0':
		HarpHost_E[m]='1e-200'
		
for n in list(HarpFungi_E.keys()):
	if HarpFungi_E[n]=='0.0':
		HarpFungi_E[n]='1e-200'

for i in list(Harp_Host.keys()):
	if i not in list(Harp_Fungi.keys()):
		abHGT.append(i)
	
	else:
		if math.fabs(math.log(float(HarpHost_E[i]))) > 2 * math.fabs(math.log(float(HarpFungi_E[i]))):
			fineHGT.append(i)
		else:
			pass

x=open('CoarseFilted_HGTgene', 'a')
y=open('FineFilted_HGTgene', 'a')
z=open('Absolute_HGTgene', 'a')
x.write('\n'.join(coarseHGT[0:]) + '\n')
x.close()
y.write('\n'.join(fineHGT[0:]) + '\n' + '\n'.join(abHGT[0:]) + '\n')
y.close()
z.write('\n'.join(abHGT[0:]) + '\n')
z.close();


