# -*- coding: utf-8 -*-
# /!\ Warning : ONLY WORKS IN Python 2.7

#To run this script it is necessary to install Heliquest
#Heliquest is available in standalone version on request to the author. R. Gautier March 2017 (https://heliquest.ipmc.cnrs.fr/)
#Please modify the Heliquest path 

HELIQUEST_PATH='/data/garrido/Logiciels/Heliquest'


import os
import pandas as pd


##########################################################################################
#SCRIPT for compile the result of the script New_All_Heliquest and find for each sequence 
#the sub sequence of lenght LEN with the best hydric moment and the best hydrophobe moment 
#Take files In_File/ResultHeliquest
#Return files in New3Helix


LEN=9 #Lenght of subsequence


def Heliquest(sequence):	
	"""
	Perform Heliquest analysis on sequence enter and return the amphipatic moment, the hydrophobic and hydrophile face of the predicted helix
	"""
	fichier = open('temp'+sequence+'.fasta', "w")
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

def All_n_seq(sequence,n):
	"""	
	From a sequence returns two lists
	- Seqs all subsequence with lenght between minimal lenght and lenght of sequence 
	- DL tuple list with start and length of all sub sequences
	"""
	Seqs = []
	D = []
	for j in range(len(sequence)-n+1):
		Seqs.append(sequence[j:j+n])
		D.append(j)
	return Seqs,D

def Best_2Helix(sequence,n):
	"""
	Find the sub sequence of lenght n in sequence wich have the best 
	hydrophobic moment and the one wich have the best hydrophobicity
	"""
	Seqs,DL=All_n_seq(sequence,n)
	MH=-9
	H=-9
	for i in range(len(Seqs)):
		 MomtHyd,FaceHydro,Hydro=Heliquest(Seqs[i])
		 if MomtHyd > MH :
		 	MH=MomtHyd
		 	dmh=i
		 if Hydro > H :
		 	H=Hydro
		 	dh=i
	return(MH,H,dmh,dh)

path_heliquest = ""

out_path = '../../New3Helix/'

for file in os.listdir(path_heliquest): 
	if file[-3:] == 'csv' and file[3]!='_':
		BD=pd.read_csv(path_heliquest+file,sep='\t')

		BDres=pd.DataFrame(columns=['PaaHheliquest','Paa.pre.Hheliquest','Paa.post.Hheliquest','Hnetcharge','Seqnetcharge','Hydro','Paa.pre.Hhydro','Paa.post.Hhydro','HydrMom','Paa.pre.Hhydrmom','Paa.post.Hhydrmom'])

		for i in range(len(BD)):
			print(str(i)+'/'+str(len(BD)))
			seq=BD.at[i,'sequence']
			if BD.at[i,'Face hydrophobe']!='[]':
				paaheli=sum(eval(BD.at[i,'Longueur']))/float(len(seq))
				paapreh=eval(BD.at[i,'Debut'])[0]/float(len(seq))
				paaposth=(len(seq)-(eval(BD.at[i,'Fin'])[-1]+1))/float(len(seq))
				print(paapreh)
			else:
				paaheli=0
				paapreh=0
				paaposth=0
			result=Best_2Helix(seq,LEN)
			hydro=result[1]
			paaprehyd=result[3]/float(len(seq))
			paaposthyd=(len(seq)-result[3]-12)/float(len(seq))
			hydmom=result[0]
			paaprehm=result[2]/float(len(seq))
			paaposthm=(len(seq)-result[2]-12)/float(len(seq))
			Hnetcharge=BD.at[i,'Helice'].count('R')+BD.at[i,'Helice'].count('K')-BD.at[i,'Helice'].count('D')-BD.at[i,'Helice'].count('E')
			Seqnetcharge=BD.at[i,'sequence'].count('R')+BD.at[i,'sequence'].count('K')-BD.at[i,'sequence'].count('D')-BD.at[i,'sequence'].count('E')
			BDres.loc[i]=[paaheli,paapreh,paaposth,Hnetcharge,Seqnetcharge,hydro,paaprehyd,paaposthyd,hydmom,paaprehm,paaposthm]
		print(out_path+'NewThreeHelix_'+file[15:])
		BDres.to_csv(out_path+'NewThreeHelix_'+file[15:],sep=',',index=False)
		print('Nombre dhelice : ', len(BD[BD['Helice']!='[]']))

