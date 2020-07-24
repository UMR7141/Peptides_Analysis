# -*- coding: utf-8 -*-
from math import *
import numpy as np
import pandas as pd
import os
import io
import sys
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from sklearn import svm
from sklearn import manifold
from sklearn import cluster, svm, metrics
from sklearn.linear_model import ElasticNet
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import classification_report,silhouette_score,silhouette_samples
from sklearn.model_selection import cross_val_score,KFold
from sklearn.cluster import KMeans,AgglomerativeClustering
from sklearn.feature_selection import VarianceThreshold,SelectFromModel, SelectFdr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet,fcluster
import scipy.stats as stats
from scipy.spatial import distance, ConvexHull
from colour import Color
from skbio import DistanceMatrix
from skbio.tree import nj




##############################
# Initialization of data set #
##############################

def Init_DFs(Ps,Nter=0):
	"""
	Initialize Data from the files and parametres given in data.py

	IN : 
		Ps -> List of desired data (initialized in data.py)
		Nter -> 0 : use all sequences / 15 : use only 15 first aa / -15 : use only 15 last aa

	OUT : 
		BDs -> Lists of Dataframes containing the names of the peptides and sequences
		couleur -> List of colors associated with each grops of data
		Label -> List of labels associated with each group
		Seq -> List of column names containing the sequences associated with each dataframe of each group
		Name -> Lists of column names containing the names of the peptides associated with each dataframe of each group
		Each of the above lists are ordered in the same way (same as Ps list)
	"""
	BDs = [] 
	couleur = []
	Label = []
	Seq = [] 
	Name = [] 
	for P in Ps : 
		print(P)
		BDs.append(pd.read_csv(P[0],sep=P[1]))
		couleur.append(P[2])		
		Label.append(P[3])
		Seq.append(P[4])
		Name.append(P[5])
	return BDs,couleur,Label,Seq,Name

def Load_DFs(Ps):	
	"""
	Initialize Data from the files saved last time and parametres given in data.py

	IN : 
		Ps -> List of desired data (initialized in data.py)

	OUT : 
		BDs -> Lists of Dataframes containing the names of the peptides and sequences
		couleur -> List of colors associated with each grops of data
		Label -> List of labels associated with each group
		Seq -> List of column names containing the sequences associated with each dataframe of each group
		Name -> Lists of column names containing the names of the peptides associated with each dataframe of each group
	"""
	couleur = []
	Label = []
	Seq = []
	Name = []
	for P in Ps :
		couleur.append(P[2])		
		Label.append(P[3])
		Seq.append(P[4])
		Name.append(P[5])
	BDs = All_New(Label) #Reload data with file save in New_csv 
	return BDs,couleur,Label,Seq,Name


##################
# Deleting files #
##################

def Clean_director():
	"""
	Delete all files from New_csv, Z_scale_Wold and ACC folders 
	"""
	for f in os.listdir('New_csv'):
		os.remove('New_csv/'+f)
	for f in [x for x in os.listdir('ACC/') if str(x[-3:])=='csv']:
		os.remove('ACC/'+f)
	for f in os.listdir('Z_scale_Wold'):
		os.remove('Z_scale_Wold/'+f)


#############
# CSV files #
#############

def Crea_New(BDs,Label):
	"""
	Create csv backup files of sequences in New_csv folder. 
	
	IN :
		BDs -> Lists of Dataframes containing the names of the peptides and sequences
		Label -> List of labels associated with each group

	"""
	for i in range(len(BDs)):
		BDs[i].to_csv('New_csv/'+Label[i]+'.csv',index=False,sep='\t')
 
def All_New(Label):
	"""
	Read file in New_csv folder
	(use in function Load_DFs

	IN : 
		Label -> List of labels associated with each group
	OUT :
		BDs -> Lists of Dataframes containing the names of the peptides and sequences
	"""
	BDs=[]
	for i in range(len(Label)):
		bd=pd.read_csv('New_csv/'+Label[i]+'.csv',sep='\t')
		BDs.append(bd)
	return BDs


###################################
# Z-scale factor Hellbergerg Wold #
###################################

def Run_Z(Label,Seq) :
	"""
	Creation of the csv files containing the 3 average z-scale values for each peptides of each of the groups in
	the Z_scale_Wold folder
	"""
	for i in range(len(Label)):
		os.system('Rscript --vanilla Z_scale_Wold.R -f '+'New_csv/'+Label[i]+'.csv -a '+Seq[i]+' -o Z_scale_Wold/z_'+Label[i]+'.csv')

def All_Z(Label):
	"""
	Reads the csv files from the Z_scale_Wold folder of interest groups and returns the Z-scale Dataframes list
	"""
	Z=[]
	for i in range(len(Label)):
		z=pd.read_csv('Z_scale_Wold/z_'+Label[i]+'.csv',sep=',')
		Z.append(z)
	return Z


##################
# ACC on Z-scale #
##################

def Run_ACC(BDs,Label,Seq,lag) :
	"""
	Creation of csv files for each group containing the lag*9 termes ACC z-scale for each peptide  
	"""
	for i in range(len(BDs)):
		os.system('Rscript --vanilla acc.R -f '+'New_csv/'+Label[i]+'.csv -a '+Seq[i]+' -l '+str(lag)+' -o ACC/acc_'+Label[i]+'.csv')

def All_ACC(BDs,Label):	
	"""
	Read csv files for each group containing the lag*9 termes ACC z-scale for each peptide  
	"""
	ACC=[]
	for i in range(len(BDs)):
		acc = pd.read_csv('ACC/acc_'+Label[i]+'.csv',sep=',')
		ACC.append(acc)
	return ACC


######################
# Update settings.py #
######################

def MaJ_Pancien_Param():
	"""
	Updating P_ancien of the parametres.py file
	P_ancien indicates which data was used the last time the program was used
	and therefore indicates what data is present in the New_csv and ACC folders
	"""
	don = open('data.py','r')
	param = open('settings.py','r')
	param = param.readlines()
	new_param = open('settings.py','w')
	for l in param:
		if 'P_ancien' in l:
			for ll in don :
				if ll[0:2]=='Ps':
					new_param.write('P_ancien = '+ll.replace(' ','')[3:])
		else :
			new_param.write(l)
	don.close()
	new_param.close()

def MaJ_Bornes_Param(Inf,Sup):
	"""
	Update of the Inf and Sup parameters of the parametres.py file
	Inf and Sup correspond to the minimum and maximum length for a peptide in order to be taken into account in the analysis
	"""
	param = open('settings.py','r')
	param = param.readlines()
	new_param = open('settings.py','w')
	for l in param:
		if 'Inf' in l:
			new_param.write('Inf = '+str(Inf)+'\n')
		elif 'Sup' in l:
			new_param.write('Sup = '+str(Sup)+'\n')
		else :
			new_param.write(l)
	new_param.close()

def MaJ_lag_Param(lag):
	"""
	Update lag settings of the parametres.py file
	lag is the number of neighbors taken into account for the calculation of auto-cross covariance (ACC)
	"""
	param = open('settings.py','r')
	param = param.readlines()
	new_param = open('settings.py','w')
	for l in param:
		if 'lag' in l:
			new_param.write('lag = '+str(lag)+'\n')
		else :
			new_param.write(l)
	new_param.close()


########################
# Displays data groups #
########################

def Aff_group(Label):
	"""
	Displays in the terminal the index associated with each group present in the analysis
	"""
	for i in range(len(Label)):
		print (i, ' : ',Label[i])


####################################
# Selection of sequences by lenght #
####################################

def Selec(BD,seq,inf,sup):
	"""
	Modify a dataframe to remove peptides that are too long or too short
	Using the Inf and Sup parameters.
	Returns the updated dataframes
	"""
	I=[]
	for i in range(len(BD)):
		if len(str(BD.at[i,seq])) >= inf and len(str(BD.at[i,seq])) <= sup :
			I.append(i)
	BD = BD.loc[I,:]
	BD.index = range(len(BD))
	return BD

def All_Selec(inf,sup,BDs,Seq):	
	"""
	Modify the dataframe in the BD list to remove peptides that are too long or short
	compared to the Inf and Sup parameters.
	Returns the list of updated dataframes
	"""
	for i in range(len(BDs)):
		BDs[i] = Selec(BDs[i],Seq[i],inf,sup)
	return BDs


###############################
# Display number of sequences #
###############################

def Nb_seq(BDs,Label):
	"""
	Display on terminal the number of peptides in each groups 
	"""
	for i in range(len(BDs)):
		print (Label[i],' : ',len(BDs[i]))


################################
# Sequence length distribution #
################################

def Add_data(BD,Seq,Label):
	"""
	Create a data frame containing three columns: label, sequences, and length
	From DB a dataframe, Seq the name of the column of this data frame containing
	the sequences and Label the label that we want for these sequences.
	"""
	S=[]
	L=[]
	for i in range(len(BD)):
		S.append(BD.at[i,Seq])
		L.append(len(BD.at[i,Seq]))
	return pd.DataFrame({'Label' : [Label]*len(BD),'Sequence' : S,'Len' : L})

def palette_couleur(Indices,Label,couleur):
	"""
	Create a color dictionary associated with their label according to the list of Indices given
	"""
	pal={}
	for i in Indices:
		pal[Label[i]]=couleur[i]
	return pal 

def violinplot_len(Indices,BDs,Seq,Label,couleur):
	"""
	Graphical display of length distribution for each group under violinplot
	"""
	DATA = pd.concat([Add_data(BDs[i],Seq[i],Label[i]) for i in Indices])
	sns.violinplot(data=DATA,y = 'Label',x = 'Len',orient="h",cut=0,palette=palette_couleur(Indices,Label,couleur))#,inner='stick') 
	plt.savefig('Out_File/Out_Length/ViolinPlot_SeqLenght.svg')
	plt.close()
	print('The violinplot of lenght distribution has been recorded: Out_File/Out_Lenght/ViolinPlot_SeqLenght.svg')





#############################
# AA Frequency distribution #
#############################

AA="GAVLMIPFYWSTCNQDEKRH" #List of 20 amino acid

def fre_aa_seq(seq):
	"""
	Calculates the frequency of each of the 20 aa in a given peptide sequence
	Return a list of the 20 frequency ordoned as the AA list
	"""
	F=[0]*20
	for a in seq:
		try :
			F[AA.index(a)]+=1
		except :
			print('Warning the sequence contains an unrecognized aa : ',a)
	return [x/float(sum(F)) for x in F]


def L_fre_aa(BD,Seq):
	"""
	In : BD -> Dataframe of peptide
		 Seq -> Index in dataframe of peptide sequence

	OUT : LF -> List of all list of aa frequency one list for each peptide in BD
	"""
	LF=[]
	for i in range(len(BD)):
		LF.append(fre_aa_seq(BD.at[i,Seq]))
	return LF


def test_kruskal(LFs,Indices,Label,couleur): 
	"""
	
	"""
	dfmi = pd.DataFrame(columns=pd.MultiIndex.from_product([[x for x in AA],[Label[x] for x in Indices]]),index=[Label[x] for x in Indices])
	
	for a in range(len(AA)):
		for g1 in Indices :
			for g2 in Indices[g1:]:
				if max([LFs[Indices.index(g1)][i][a] for i in range(len(LFs[Indices.index(g1)]))])!=0 or max([LFs[Indices.index(g2)][i][a] for i in range(len(LFs[Indices.index(g2)]))])!=0:# pour contrer l'erreur : ValueError: All numbers are identical in kruskal car toutes les valeur sont 0 
					pk=stats.kruskal([LFs[Indices.index(g1)][i][a] for i in range(len(LFs[Indices.index(g1)]))],[LFs[Indices.index(g2)][i][a] for i in range(len(LFs[Indices.index(g2)]))])[1]
					dfmi.at[Label[g1],(AA[a],Label[g2])]= pk	
	dfmi.to_csv('Out_File/Out_Freq/Out_File_Test_kruskal.csv',sep='\t',index=True)
	print('The file of multiple kruskal statistical test has been recorded: Out_File/Out_Freq/Out_File_Test_kruskal.csv')
	A=[]
	F=[]
	L=[]
	for i in range(len(LFs)):
		for j in range(len(LFs[i])):
			for a in range(len(LFs[i][j])):
				L.append(Label[Indices[i]])
				F.append(LFs[i][j][a])
				A.append(AA[a])
	DATA = pd.DataFrame({'AA':A,'freq':F,'Label':L})
	plt.figure(1, figsize=(5*len(Indices), 20))
	sns.boxplot(data = DATA , x = 'AA', y = 'freq', hue='Label', fliersize = 0,palette=palette_couleur(Indices,Label,couleur))
	plt.savefig('Out_File/Out_Freq/HistogramAAfrequency.svg')
	plt.close()
	print('The figure of histogram frequency has been recorded: Out_File/Out_Freq/HistogramAAfrequency.svg\n\n')
	

def plot_dist_feq_aa (Indices,BDs,Label,Seq,couleur):
	"""
	plot the aa frequency distribution of all the groups in the list Indices
	"""
	e = 0
	LFs=[]
	for i in Indices :
		LF =  L_fre_aa(BDs[i],Seq[i])

		LFs.append(LF)
		Mean_LF = [np.mean([LF[j][x] for j in range(len(LF))]) for x in range(len(LF[0]))]
		Std_LF = [np.std([LF[j][x] for j in range(len(LF))]) for x in range(len(LF[0]))]
		e += 1/float(len(Indices)+1)
	test_kruskal(LFs,Indices,Label,couleur)




##################################
# Topology heat map distribution #
##################################


def Dist_Heat(Indice,m,BDs,Seq,Label):
	"""
	Makes a heat map of the distribution of aa according to their relative position on the peptide in a specific group indicate in Indice
	Result figure : Out_File/HeatMap_AA_pos.svg
	"""
	Mat = np.zeros((20,1100))
	for i in range(len(BDs[Indice])):
		seq = BDs[Indice].at[i,Seq[Indice]]
		for a in range(len(seq)) :
			for j in range (int(float(a-0.5)/(len(seq)-1)*999),int(float(a+0.5)/(len(seq)-1)*999)):
				if j <1000 and j>0:
					try :
						Mat[AA.index(seq[a])][j]+=1
					except :
						print('Attention la sequence contient un aa non reconnu : ',seq[a])
	Mat=Mat/np.sum(Mat,axis=0)
	for j in range(20):
		M = np.sum([x for x in Mat[j][:1000] if str(x)!='nan'])/1000
		for i in range(1050,1100):
			Mat[j][i]=M
	Mat = pd.DataFrame(Mat)
	plt.figure(1, figsize=(25, 20))
	ax = sns.heatmap(Mat,cmap="Blues",vmax=m)
	ax.yaxis.set_ticklabels("GAVLMIPFYWSTCNQDEKRH", rotation = 0, fontsize = 40, verticalalignment = 'center')
	ax.xaxis.set_visible(False)
	plt.title(Label[Indice], fontsize=35)
	plt.tick_params(axis = 'y', length = 0)
	plt.savefig('Out_File/Out_Freq/HeatMap_AA_pos_'+Label[Indice]+'.svg')
	print('The figure has been recorded: Out_File/Out_Freq/HeatMap_AA_pos_'+Label[Indice]+'.svg')
	plt.close()

def Dist_Heat_Moy(Indice,m,BDs,Seq,Label):
	"""
	Make a heat map of the mean frequency of each aa observed in serval group indicate in Indice list 
	Result figure : Out_File/HeatMap_AA_mean.svg
	"""
	Mat = np.zeros((20,len(Indice)))
	for ind in Indice :
		for i in range(len(BDs[ind])):
			seq = BDs[ind].at[i,Seq[ind]]
			for a in range(len(seq)) :
				try :
					Mat[AA.index(seq[a])][Indice.index(ind)]+=1
				except :
					print('Ya un ',a,' dans',Label[ind])
				

	Mat=Mat/np.sum(Mat,axis=0)

	Mat = pd.DataFrame(Mat,index=['G','A','V','L','M','I','P','F','Y','W','S','T','C','N','Q','D','E','K','R','H'],columns=[Label[x] for x in Indice])
	plt.figure(1, figsize=(5*len(Indice), 20))
	ax = sns.heatmap(Mat,cmap="Blues",vmax=m)
	ax.yaxis.set_ticklabels("GAVLMIPFYWSTCNQDEKRH", rotation = 0, fontsize = 40, verticalalignment = 'center')
	ax.xaxis.set_ticklabels([Label[x] for x in Indice], rotation = 0, fontsize = 35, verticalalignment = 'center')

	plt.tick_params(axis = 'y', length = 0)
	plt.savefig('Out_File/Out_Freq/HeatMap_AA_pos.svg')
	print('The figure has been recorded: Out_File/Out_Freq/HeatMap_AA_pos.svg')
	plt.close()








###################################
# Z-scale factor Hellbergerg Wold #
###################################

def plot_Zscale_1D(Indices,Label,couleur,Z):
	"""
	Plot mean of z-scale 
	"""
	ax=plt.gca()
	plt.figure(1, figsize=(40, 40))
	for j,i in enumerate(Indices):
		plt.plot(Z[i]['V1'],[1+(0.1*j)]*len(Z[i]['V1']),'o',color=couleur[i])
		plt.plot(Z[i]['V2'],[2+(0.1*j)]*len(Z[i]['V2']),'o',color=couleur[i])
		plt.plot(Z[i]['V3'],[3+(0.1*j)]*len(Z[i]['V3']),'o',color=couleur[i])
	ax.set_yticks([1,2,3])
	ax.yaxis.set_ticklabels(['Z1','Z2','Z3'], rotation = 0, fontsize = 25, verticalalignment = 'center')
	Legend = [mpatches.Patch(color=couleur[i], label=Label[i]) for i in Indices]
	plt.legend(handles=Legend)
	plt.savefig('Out_File/Out_Zscale/Zscale1D.svg')
	plt.close()
	print('The figure has been recorded: Out_File/Out_Zscale/Zscale1D.svg')

def plot_Zscale(Indices,Label,couleur,Z):
	"""
	Plot Mean Z-scale in 3 dimensions
	"""
	fig = plt.figure(1, figsize=(4, 3))
	plt.clf()
	ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

	for i in Indices:
		ax.scatter(Z[i]['V1'], Z[i]['V2'], Z[i]['V3'], cmap=plt.cm.nipy_spectral,color=couleur[i], depthshade=False)#Pour que les points change pas de couleur mettre False à depthshade
	ax.set_xlabel('z1')
	ax.set_ylabel('z2')
	ax.set_zlabel('z3')
	Legend = [mpatches.Patch(color=couleur[i], label=Label[i]) for i in Indices]
	plt.legend(handles=Legend)
	plt.show()


##########
# PCA_2D #
##########

def Corel_circle(acp,n,p,NameF):
	"""
	Create an svg file containing the PCA corel circle figure
	IN : acp -> model acp
	     n -> number of individuals
	     p -> number of feature
	     NameF -> list name of feature
	OUT : Out_File/Out_PCA/CorelCircle.svg
	"""
	fig = plt.figure(1, figsize=(8, 8))
	eigval=(n-1)/n*acp.explained_variance_
	sqrt_eigval = np.sqrt(eigval)
	corvar=np.zeros((p,len(eigval)))

	for k in range(len(eigval)):
		corvar[:,k]=acp.components_[k,:]*sqrt_eigval[k]
	for j in range(p):
		plt.scatter(corvar[j,0],corvar[j,1])
		plt.annotate(NameF[j],(corvar[j,0],corvar[j,1]))
	plt.plot([-1,1],[0,0],color='silver',linestyle='-',linewidth=1)
	plt.plot([0,0],[-1,1],color='silver',linestyle='-',linewidth=1)
	theta = np.linspace(0, 2*np.pi, 40)
	x = np.cos(theta)
	y = np.sin(theta)
	plt.plot(x, y)	
	theta = np.linspace(0, 2*np.pi, 40)
	x = 0.5*np.cos(theta)
	y = 0.5*np.sin(theta)
	plt.plot(x, y)
	print('The figure of the corel circle has been recorded: Out_File/Out_PCA/CorelCircle.svg')	
	plt.savefig('Out_File/Out_PCA/CorelCircle.svg')
	plt.close()

def Initializ_DF_Data(Indices,feature,z,l,auto):
	"""
	Initialization of the dataframe containing all vector of all peptides wanted for the PCA
	"""

	if z == 'all' : #Take all feature
		All = pd.concat([feature[i] for i in Indices])
	else : # Take only feature conserning one specified z-scale (z)
		C=[c for c in feature[Indices[0]] if 'z'+str(z)+'.lag' in c]
		All = pd.concat([feature[i][C] for i in Indices])

	if l!='all': #Keep only feature conserning one specified lag (l)
		for x in All:
			if 'lag'+str(l) not in x:
				del All[x]

	if auto==1: #Delete feature of Cross-covariance 
		for x in All:
			print(x,x.count('.'))
			if x.count('.') != 1 :
				del All[x]
	return (All)

def PCA_2D(Indmod,Indices, feature,BDs,Label,couleur,Name,Seq,Nam,air,z='all',l='all',L_Gros=0,auto=0,InArea=0):
	"""
	Main fonction for process PCA analyse
	IN : 
		Indmod -> Index of group for process the pca model
		Indices -> Index of group for visualization
		feature -> List of dataframe of each group containing all vector utilized for PCA
		c_coul -> List of color of each group
		BDs -> Lisr of dataframe of each group containing sequence and name
		Label -> List of Label of each group
		couleur -> List of color of each group
		Name -> List of index of dataframe in BDs corresponding to colums containing name of peptide
		Seq -> List of index of dataframe in BDs corresponding to colums containing sequence of peptide
		Nam -> If Nam = 1 Display Name of peptide on PCA figure PCA.svg
		air -> If air = 1 Display convexe area of each group of peptide on PCA figure PCA.svg
		z -> If z='all' use all ACC term / If z=1 or 2 or 3 use only ACC term concernig one specific z-scale
		l -> If l='all' use all ACC term / If l=n   use only ACC term concernig one specific lag
		L_Gros -> If =0 all peptides are represented with same point else precise index of group in list for plot with bigger point 
		auto -> If =0 use all ACC term / If =1 use only ACC terme of auto cross covariance not cross covariance
	OUT : Files in Out_File/Out_PCA/
	"""

	#=======INITIALIZATION=======

	if Indmod != 0: # To make the ACP model on a group other than the one represented
		# Initialization of matrix of all learned vectors to fit the pca
		Allmod = Initializ_DF_Data(Indmod,feature,z,l,auto)
		
	# Initialization of matrix of all represented vectors
	All = Initializ_DF_Data(Indices,feature,z,l,auto)

	Name_feature=All.columns

	# Initialization of class vector 
	Y=[]
	for i in range(len(Indices)):
		Y = Y + [i]*len(feature[Indices[i]])


	if Indmod == 0: # Same Lerned and represented ?
		Allmod=All



	#=======PROCESS PCA=======
	
	#Normalization
	scaler =  StandardScaler().fit(Allmod) 
	Allmod = scaler.transform(Allmod) 
	All = scaler.transform(All) 

	#Create model
	pca = PCA(n_components=2)
	pca.fit(Allmod)

	#Aply model
	X_pca = pca.transform(All)



	#=======REPRESENTATION=======

	Corel_circle(pca,All.shape[0],All.shape[1],Name_feature)
	
	fig = plt.figure(1, figsize=(20, 20))

	#Displays the variance explained by each of the axes
	print ('\n Explained variance of first component :',pca.explained_variance_ratio_[0],
		'\n Explained variance of second component :',pca.explained_variance_ratio_[1])





	if L_Gros != 0 : #To put some groups with different type points than others
		for i in range(len(Y)):
			if Indices[Y[i]] in L_Gros:#To put somme groups with bigger point
				plt.scatter(X_pca[:, 0][i], X_pca[:, 1][i],marker='o',s=60, c=couleur[Indices[Y[i]]])
			# elif Indices[Y[i]]==0:#The first groupe in Indices with cross marker
			# 	plt.scatter(X_pca[:, 0][i], X_pca[:, 1][i],marker='P',s=150, c=couleur[Indices[Y[i]]],linewidths=0.001)
			# elif  Indices[Y[i]]==1:#The second groupe in Indices with circle marker
			# 	plt.scatter(X_pca[:, 0][i], X_pca[:, 1][i],marker='o',s=70, edgecolor=couleur[Indices[Y[i]]],facecolor='w',linewidths=3)
			else :
				plt.scatter(X_pca[:, 0][i], X_pca[:, 1][i],marker='o',s=35, c=couleur[Indices[Y[i]]]) 
		
	else: 
		plt.scatter(X_pca[:, 0], X_pca[:, 1],marker='o',s=35, c=[couleur[Indices[y]] for y in Y]) 
	

	#For display name of peptides
	if Nam==1:
		ii=0
		for j in [Indices[0]]:#Displays the name of only the first group
			for i, txt in enumerate(BDs[j][Name[j]]):
				plt.annotate(str(i)+str(txt), (X_pca[:, 0][i+ii],X_pca[:, 1][i+ii]))
			ii+=i+1
	Legend = [mpatches.Patch(color=couleur[i], label=Label[i]) for i in Indices]
	plt.legend(handles=Legend,prop={'size':20})


	#For display convexe area
	if air == 1 :
		n=0
		for i in Indices:

			f=n+len(feature[i])
			points=X_pca[n:f]
			hull=ConvexHull(points)

			while len(points)>float(f-n)*0.50 and len(points)>10: #To display central area about 50% of the points
				points = np.delete(points, hull.vertices, axis=0)
				hull=ConvexHull(points)
			x=[]
			y=[]	

			for simplex in hull.simplices: # Pour tracer les contours
				plt.plot(points[simplex, 0], points[simplex, 1], couleur[i], linewidth=5)


			#For create a file with the sequences of InArea[0] and a new column that indicate if the peptide is in (1) or out (0) 
			#the convexe area of the groupe InArea[1]
			if InArea != 0 :
				if i == InArea[1] : #Index of the group to test her area
					In_Out_Convex_Area_file(BDs[InArea[0]],hull.vertices,hull.points,X_pca[:, 0],X_pca[:, 1])#Index of the group to test if their peptides are in the area

			n=f


	plt.axis(( min(min(X_pca[:, 0]),min(X_pca[:, 1]))-0.5 , max(max(X_pca[:, 0]),max(X_pca[:, 1]))+0.5 , min(min(X_pca[:, 0]),min(X_pca[:, 1]))-0.5 , max(max(X_pca[:, 0]),max(X_pca[:, 1]))+0.5))
	plt.xlabel('PC1 : '+str(pca.explained_variance_ratio_[0]))
	plt.ylabel('PC2 : '+str(pca.explained_variance_ratio_[1]))


	plt.tick_params(labelsize=16)
	plt.savefig('.svg')			
	plt.savefig('Out_File/Out_PCA/PCA.svg')
	print('The figure of the PCA has been recorded: Out_File/Out_PCA/PCA.svg')
	plt.close()


	#For create file with all the coordinate
	Label=[]
	S=[]
	for j in Indices:
		for i in range(len(BDs[j])):
			Label.append(BDs[j][Name[j]][i])
			S.append(BDs[j][Seq[j]][i])
	BD = pd.DataFrame({'Name' : Label,'Sequence' : S,'X' :X_pca[:, 0],'Y': X_pca[:, 1] })
	BD.to_csv('Out_File/Out_PCA/Coordinate.csv',index=False)
	print('The file of coordinate has been recorded: Out_File/Out_PCA/Coordinate.csv')

def In_Out_Convex_Area_file(BD,IndArea,coordArea,X,Y):
	'''
	Saves a csv file corresponding to the entry csv file with a columns in addition
	which indicates whether the peptides in the desired representation is in or out the desired convex area
	To change group peptide tested and the group of convex area need to modify the PCA_2D function 

	IndArea index lists ordered backward direction clockwise hull.vertices
	coordArea contains the coordinates in the ACP of the points in the area given by corespondance with the indices the coordinates
	'''	
	X_Area=[]
	Y_Area=[]
	for i in IndArea:
		X_Area.append(coordArea[i][0])
		Y_Area.append(coordArea[i][1])
	X_Area.append(X_Area[0])
	Y_Area.append(Y_Area[0])
	for i in range(len(BD)):
		BD.at[i,'ConvexArea']=Point_in_Convexe_Area(X_Area,Y_Area,X[i],Y[i])
	BD.to_csv('Out_File/Out_PCA/InOutConvexeArea.csv')
	print('The file of in out convexe  area has been recorded: Out_File/Out_PCA/InOutConvexeArea.csv')

def Point_in_Convexe_Area(X_Area,Y_Area,x_test,y_test):
	'''
	Determine if one point is in or out a convexe area
	X_Area : List of X coordinate of convexe area point
	Y_Area : List of Y coordinate of convexe area point
	x_test : x cordinate of tested point
	y_test : y cordinate of tested point
	'''
	SA=0 #Initialize sum of angle
	for i in range(len(X_Area)-1):
		V1=[X_Area[i]-x_test,Y_Area[i]-y_test]
		V2=[X_Area[i+1]-x_test,Y_Area[i+1]-y_test]
		SA+=degrees(acos( ((V1[0]*V2[0]) +(V1[1]*V2[1])) /  (sqrt(V1[0]**2+V1[1]**2)*sqrt(V2[0]**2+V2[1]**2)) ))
	if round(SA,0)==360:
		return(1) #If sum is 360 the point is In the area
	else :
		return(0)


	
##########
# PCA_3D #
##########

def PCA_3D(Indices, feature, BDs,Label,couleur,Name,Seq):

	#=======INITIALIZATION=======

	# Initialization of matrix of all represented vectors
	All = Initializ_DF_Data(Indices,feature,'all','all',0)

	Name_feature=All.columns

	# Initialization of class vector 
	Y=[]
	for i in range(len(Indices)):
		Y = Y + [i]*len(feature[Indices[i]])



	Allmod=All


	#=======PROCESS PCA=======
	
	#Normalization
	scaler =  StandardScaler().fit(Allmod) 
	Allmod = scaler.transform(Allmod) 
	All = scaler.transform(All) 

	#Create model
	pca = PCA(n_components=3)
	pca.fit(Allmod)

	#Aply model
	X_pca = pca.transform(All)



	#=======REPRESENTATION=======

	#Graphic representation
	fig = plt.figure(1, figsize=(4, 3))
	plt.clf()
	ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

	ax.scatter(X_pca[:, 0], X_pca[:, 1], X_pca[:, 2], c=[couleur[Indices[y]] for y in Y], cmap=plt.cm.nipy_spectral)


	Legend = [mpatches.Patch(color=couleur[i], label=Label[i]) for i in Indices]

	plt.legend(handles=Legend)	
	ax.set_xlabel('PC1 :'+str(round(pca.explained_variance_ratio_[0],3)))
	ax.set_ylabel('PC2 :'+str(round(pca.explained_variance_ratio_[1],3)))
	ax.set_zlabel('PC3 : '+str(round(pca.explained_variance_ratio_[2],3)))
	#plt.savefig('Out_File/Out_PCA/PCA3D.csv')
	plt.show()
	plt.close()
	#print('The figure of the 3D PCA has been recorded: Out_File/Out_PCA/PCA3D.csv')

	#Output file with coordinates
	Label=[]
	S=[]
	for j in Indices:
		for i in range(len(BDs[j])):
			Label.append(BDs[j][Name[j]][i])
			S.append(BDs[j][Seq[j]][i])
	BD = pd.DataFrame({'Name' : Label,'Sequence' : S,'X' :X_pca[:, 0],'Y': X_pca[:, 1],'Z': X_pca[:, 2] })
	BD.to_csv('Out_File/Out_PCA/PCA3DCoord.csv',index=False)
	print('The file of the 3D coordinate PCA has been recorded: Out_File/Out_PCA/PCA3DCoord.csv')
	




###################
# Distance Matrix #
###################



def Make_Dist_Mat(Ind,feature,BDs,Label,Name):
	'''
	Make a distance matrix 
	IN:
		Ind
		feature
		BDs
		Label
		Name
	OUT:
		Matrix distance phylip format: Out_File/MatDist
	'''	

	#Matrix distance initialization
	mat_distND=np.zeros((len(feature[Ind]),len(feature[Ind])))

	#Fill the matrix with all pair euclidean distance
	for i in range(len(feature[Ind])):
		for j in range(i+1,len(feature[Ind])):
			
			mat_distND[i][j]=distance.euclidean(feature[Ind].loc[i],feature[Ind].loc[j])
			mat_distND[j][i]=mat_distND[i][j]

	#Make Phylip file 
	f=open('Out_File/Out_Tree/MatDist.phylip','w')
	f.write('    '+str(len(mat_distND))+'\n')
	for i,l in enumerate(BDs[Ind][Name[Ind]]):

		s=''
		for x in range(len(mat_distND[i])):
			#print(mat_distND)
			s+=str(round(mat_distND[i][x],3))
			if x!= len(mat_distND[i])-1:
				s+='  '
			else:
				s+='\n'
		f.write('{:10}'.format(l[:min(len(l),9)])+s)
	f.close()
	print('The phylip file of distance matrix has been recorded: Out_File/Out_Tree/MatDist.phylip')
	return(mat_distND)

def Make_newick_Tree(Ind,feature,BDs,Label,Name):
	mat_distND=Make_Dist_Mat(Ind,feature,BDs,Label,Name)
	dm=DistanceMatrix(mat_distND,list(BDs[Ind][Name[Ind]]))

	tree = nj(dm)
	
	newick_str = nj(dm, result_constructor=str)
	f=open('Out_File/Out_Tree/Newick_tree','w')
	f.write(newick_str)
	print('The newick file of the tree has been recorded: Out_File/Out_Tree/Newick_tree')
	



###########
# K-MEANS #
###########


def Process_Kmean(Ind,feature,BDs,Label,Name,colors,k):
	'''
	Process a K-mean classification
	Ind: index of groups use for kmean
	feature: List of group vector dataframe (ex : ACC)
	BDs: list of dataframe contain Name and Sequence of all peptides
	Label: list of label
	Name: list of column of name peptide in each dataframe in BDs
	k: number of cluster 
	'''
	X=[] #Initialize list of feature vector
	Lab=[] #List of label (X order)
	Nam=[]#List of name (X order)
	for ind in Ind:
		for i in range(len(feature[ind])):
			X.append(feature[ind].loc[i])
			Lab.append(Label[ind])
			Nam.append(BDs[ind].at[i,Name[ind]])

	X = np.array(X)
	kmeans = KMeans(n_clusters=k,random_state=0,n_init=100).fit(X) 

	result = pd.DataFrame(np.zeros((len(Ind)+1,k)),columns=['Cluster '+str(x) for x in range(1,k+1)], index=[Label[x] for x in Ind]+['Silhouette coef'])


	Nam_Group=[[] for x in range(k)]
	Lab_Group=[[] for x in range(k)]
	for i in range(len(kmeans.labels_)):
		result.at[Lab[i] , 'Cluster '+str(kmeans.labels_[i]+1)]+=1
		Nam_Group[kmeans.labels_[i]].append(Nam[i])
		Lab_Group[kmeans.labels_[i]].append(Lab[i])



	
	all_silhouette = silhouette_samples(X, kmeans.labels_)
	for i in range(k):
		Sil=[all_silhouette[x] for x in range(len(all_silhouette)) if kmeans.labels_[x]==i]
		result.at['Silhouette coef','Cluster '+str(i+1)]=np.mean(Sil)
	result.at['Silhouette coef','Global']=silhouette_score(X, kmeans.labels_)




	result.to_csv('Out_File/Out_Kmean/Tab_Kmean.csv')
	AllPep=pd.DataFrame({'Name':Nam,'Label':Lab,'Cluster':kmeans.labels_+1})
	sns.catplot(data=AllPep,hue='Label',x='Cluster',kind="count",palette=[colors[x] for x in Ind])
	plt.savefig('Out_File/Out_Kmean/ByCluster.svg')
	plt.close()
	sns.catplot(data=AllPep,hue='Cluster',x='Label',kind="count")
	plt.savefig('Out_File/Out_Kmean/ByClass.svg')
	plt.close()
	AllPep.to_csv('Out_File/Out_Kmean/All_Pep_Cluster.csv',index=False)
	os.system('Rscript MultChi2.r')

