# -*- coding: utf-8 -*-
# /!\ Warning : ONLY WORKS IN Python 2.7

#To run this script it is necessary to install Heliquest
#Heliquest is available in standalone version on request to the author. R. Gautier March 2017 (https://heliquest.ipmc.cnrs.fr/)
#Please modify the Heliquest path 

HELIQUEST_PATH='/data/garrido/Logiciels/Heliquest'


import pandas as pd
import os
from math import *
OY={'A':0,'L':0,'I':0,'V':0,'M':0,'P':0,'F':0,'W':0,'Y':0,'E':1,'D':1,'K':1,'R':1,'S':1,'T':1,'N':1,'Q':1,'H':1,'G':2,'C':3}#Dico 0:Hydrophobe 1:Hydrophyle 3:exception
Face = {0:9,9:0,1:10,10:1,2:11,11:2,3:12,12:3,4:13,13:4,5:14,14:5,6:15,15:6,7:16,16:7,8:17,17:8}


LEN=9 #Minimum helix length

def All_seq(sequence):
	"""
	From a sequence returns two lists
	- Seqs all subsequence with lenght between minimal lenght and lenght of sequence 
	- DL tuple list with start and length of all sub sequences
	"""
	Seqs = []
	DL=[]
	for i in range(LEN,len(sequence)+1)[::-1]:
		for j in range(len(sequence)-i+1):
			Seqs.append(sequence[j:j+i])
			DL.append((j,i))
	return Seqs,DL




def Heliquest(sequence):
	"""
	Perform Heliquest analysis on sequence enter and return the amphipatic moment, the hydrophobic and hydrophile face of the predicted helix
	"""
	fichier = open("temp"+sequence+".fasta", "w")
	fichier.write(">Seq\n"+sequence)
	fichier.close()	
	os.system('python '+HELIQUEST_PATH+'/HeliQuest_analysis.py -i temp'+sequence+'.fasta  -o res_temp'+sequence+' -t full')
	Res_hel = pd.read_csv('res_temp'+sequence+'.csv',sep=';')
	MomtHyd = Res_hel.at[0,'  Hyd.Moment ']
	FaceHydro = Res_hel.at[0,' Hydrophobic Face']
	Hydro = Res_hel.at[0,' Hydrophobicity']	
	os.remove('temp'+sequence+'.fasta')
	os.remove('res_temp'+sequence+'.csv')
	return(MomtHyd,FaceHydro,Hydro)


def Helical_voisin(sequence):
	'''
	Takes as argument the aa sequence of the helix and returns
	a list of 18 character strings corresponding to the 18 positions 
	different from the helix
	'''
	Order = [0,11,4,15,8,1,12,5,16,9,2,13,6,17,10,3,14,7]
	Seq_order = ['']*18
	for i in range(len(sequence)):
		Seq_order[Order.index(i%18)]+=sequence[i]
	return(Seq_order)

#print(Helical_voisin(sequence))

def Find_best_HF(sequence):
	'''
	Find the best possible hydrophobic face on a sequence
	Return these hydrophobic face and the associated hydrophyle face 

	'''
	Seq_order=Helical_voisin(sequence)
	Seq_order_OY=[] #Same List of 18 lists as Seq_order but each aa is replaced by a number corresponding to hydrophobic hydrophilic .
	for pos in Seq_order:
		Seq_order_OY.append([OY[x] for x in pos])
	HF=''
	i=0
	while i < 18 : #All possible starts 
		deb=i
		hf=''
		if 1 in Seq_order_OY[i] or 3 in Seq_order_OY[i] or Seq_order_OY[i]==[]: 
			i+=1
		else :
			j=0
			while j < 18 and 1 not in Seq_order_OY[(i+j)%18] and 3 not in Seq_order_OY[(i+j)%18]:
				longueur=j+1
				hf+=Seq_order[(i+j)%18]
				j+=1

			#Do not take the empty positions at the end of the face.
			while Seq_order_OY[(i+longueur)%18]==[]:
				longueur-=1

			#Special case of wisteria if it is at the end of the face, it is removed from the face.
			while 2 in Seq_order_OY[deb] :
				hf=hf.replace('G','',Seq_order[i].count('G'))  # Pour enlever le ou les G si il font partie de la premiere position mais garde les autres
				if sum(Seq_order_OY[i])/len(Seq_order_OY[i])==2: #Regarde si yavait que des G et dans ce cas diminue le debut et sinon peut partir
					deb=(deb+1)%18
				else:
					break
			while 2 in Seq_order_OY[(deb+longueur-1)%18]:
				hf=hf[::-1].replace('G','',Seq_order[i].count('G'))[::-1]  #To remove the G if they are part of the first position but keep the others.
				if sum(Seq_order_OY[i])/len(Seq_order_OY[i])==2:
					longueur-=1
				else:
					break


			if len(hf)>len(HF):
				HF=hf
				Deb=deb
				Longueur=longueur
			i+=j
	if HF!='':
		face_hydrophyle=Find_YF(Seq_order,Deb,Longueur)
		if Longueur<3 or len(HF)<5 or face_hydrophyle==0 : #If the hydrophobic side is not wide enough (less than three positions / 18 ) or if the hydrophilic side does not contain enough polar aa 
			return(0)
	else:
		return(0)
	return (HF,face_hydrophyle)





def Find_YF(Seq_order,deb,longueur):
	'''
	Check that the residues in front of the hydrophobic side are polar.
	Takes a window of three or four residues in front depending on whether 
	the length of the hydrophobic face is even or odd
	'''
	deb = Face[deb]

	#The hydrophilic side : all residues in front of the hydrophobic side
	face_hydrophyle = ''
	for i in range(deb,deb+longueur):
		face_hydrophyle += Seq_order[i%18] 

	#Reduce the size of the hydrophilic window to verify that the residues in the center (3 or 4) are all hydrophilic.
	while longueur>4:
		deb+=1
		longueur-=2

	#Check that the residues in the centre (3 or 4) are all hydrophilic.
	for i in range(deb,deb+longueur):
		for a in Seq_order[i%18]:
			if a not in 'AEDKRSTNQHG': 
				return (0)
	return (face_hydrophyle)


def Best_Helix(sequence):
	'''
	Return
	The sequence of predicted helix, the hydrophobic face, the hydrophyle face, start of helix, end of helix, lenght of helix, moment hydrique of helix, hydrophobicity of helix
	'''
	print(sequence)
	Seqs,DL = All_seq(sequence)
	face_hydrophobe='NONE'
	i=0
	while (face_hydrophobe=='NONE') and i<len(Seqs):

		Result=Find_best_HF(Seqs[i])
		if Result!=0:
			face_hydrophobe=Result[0]
		i+=1
	if i!=len(Seqs)-1:
		i-=1
	if 	Result!=0:
		debut = DL[i][0]
		longueur=DL[i][1]
		moment_hydrique,FaceHydro,Hydro=Heliquest(Seqs[i])
		fin=debut+longueur-1
		return (Seqs[i],Result[0],Result[1],debut,fin,longueur,moment_hydrique,Hydro)
	else:

		return(0)

#print(Find_best_HF(sequence))
#print(Best_Helix(sequence))



def Result_Heliquest(BD,Seq,Out):
	'''
	BD : dataframe with sequence to predict, 
	Seq : column name of sequences
	'''
	D=[]#Start helix
	F=[]#End helix
	S=[]#Complete peptide sequence
	H=[]#Selected helix
	MH=[]#Amphipathic moment
	Hy=[]#Hydrophobicity
	FHO=[]#Face hydrophobe
	FHY=[]#Face hydrophyle
	L=[]#Lenght helix
	for i in range(len(BD)):
		print(str(float(i)/float(len(BD))*100.)+'% ')
		Sequence=BD.at[i,Seq]
		d=[]
		f=[]
		h=[]
		mh=[]
		hy=[]
		fho=[]
		fhy=[]
		l=[]
		seqs=[]#Liste des bout de sequence a analyser
		debs=[]#Liste des indices du début des seq a analyser as rapport à la séquence complete
		if len(Sequence)>LEN-1:
			seqs.append(Sequence)
			debs.append(0)
		while seqs!=[]:
			#for i in range(len(seqs)):
				i=0  # i is always 0 because we clear the list of sequence seqs to check as we go along 
				Result = Best_Helix(seqs[i])
				if Result!=0:
					d.append(Result[3]+debs[i])
					f.append(Result[4]+debs[i])
					h.append(Result[0])
					mh.append(Result[6])
					hy.append(Result[7])
					fho.append(Result[1])
					fhy.append(Result[2])
					l.append(Result[5])

					if Result[3]>LEN-1:
						debs.append(debs[i])
						seqs.append(seqs[i][:int(Result[3])])

					if len(seqs[i])-Result[4]-1>LEN-1:
						debs.append(f[-1]+1)
						seqs.append(seqs[i][Result[4]+1:])

				del(seqs[i])
				del(debs[i])

		S.append(Sequence)
		D.append(d)
		F.append(f)
		H.append(h)
		MH.append(mh)
		Hy.append(hy)
		FHO.append(fho)
		FHY.append(fhy)
		L.append(l)
	pd.DataFrame({'Name' : range(len(BD)), 'sequence' : S, 'Debut' :D, 'Fin' : F, 'Hydrophobicity' : Hy,
 	 'Hyd.Moment' : MH, 'Face hydrophobe' : FHO, 'Face hydrophyle' : FHY,'Helice' : H, 'Longueur' : L}).to_csv(Out,index=False,sep='\t')


# Result_Heliquest(pd.read_csv('../HA-RAMP/Group1.csv',sep='\t'),'Sequence','result.heliquest.Group1.csv')
# Result_Heliquest(pd.read_csv('../HA-RAMP/Group2.csv',sep='\t'),'Sequence','result.heliquest.Group2.csv')
# Result_Heliquest(pd.read_csv('../TP/mTP.csv',sep='\t'),'Peptide','result.heliquest.mTP.csv')
# Result_Heliquest(pd.read_csv('../SP/SP.csv',sep='\t'),'Sequence','result.heliquest.SP.csv')
# Result_Heliquest(pd.read_csv('../TP/cTP.csv',sep='\t'),'Peptide','result.heliquest.cTP.csv')
# Result_Heliquest(pd.read_csv('../Other/Random.csv',sep=','),'Sequence','result.heliquest.Random.csv')
# Result_Heliquest(pd.read_csv('../TP/mTPChlamydomonas.csv',sep='\t'),'Peptide','result.heliquest.mTPChlamydomonas.csv')
# Result_Heliquest(pd.read_csv('../TP/mTPArabidopsis.csv',sep='\t'),'Peptide','result.heliquest.mTPArabidopsis.csv')
# Result_Heliquest(pd.read_csv('../TP/mTPHomo.csv',sep='\t'),'Peptide','result.heliquest.mTPHomo.csv')
# Result_Heliquest(pd.read_csv('../TP/mTPCerevisae.csv',sep='\t'),'Peptide','result.heliquest.mTPCerevisae.csv')
# Result_Heliquest(pd.read_csv('../SP/bSP.csv',sep='\t'),'Sequence','result.heliquest.bSP.csv')
# Result_Heliquest(pd.read_csv('../SP/tSP.csv',sep='\t'),'Sequence','result.heliquest.tSP.csv')
# Result_Heliquest(pd.read_csv('../SP/eSP.csv',sep='\t'),'Sequence','result.heliquest.eSP.csv')
# Result_Heliquest(pd.read_csv('../TP/cTPArabidopsis.csv',sep='\t'),'Peptide','result.heliquest.cTPArabidopsis.csv')
# Result_Heliquest(pd.read_csv('../TP/cTPChlamydomonas.csv',sep='\t'),'Peptide','result.heliquest.cTPChlamydomonas.csv')
# Result_Heliquest(pd.read_csv('../Other/Cyclotide.csv',sep='\t'),'Sequence','result.heliquest.Cyclotide.csv')