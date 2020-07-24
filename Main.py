# -*- coding: utf-8 -*-

from data import *
from settings import *
from functions import *




#########################
# Initialization of DFs #
#########################

if P_ancien != Ps :# Check if data change, if is it update all files 
	Clean_director()
	BDs,couleur,Label,Seq,Name = Init_DFs(Ps)
	BDs=All_Selec(Inf,Sup,BDs,Seq)
	Crea_New(BDs,Label)
	Run_ACC(BDs,Label,Seq,lag)
	ACC=All_ACC(BDs,Label)
	Run_Z(Label,Seq)
	Z=All_Z(Label)
	MaJ_Pancien_Param()

else :# Don't change reload data with saved files
	BDs,couleur,Label,Seq,Name = Load_DFs(Ps)
	ACC=All_ACC(BDs,Label)
	Z=All_Z(Label)



########
# Main #
########

A = 1
while A != 0:

	print('\n  MENU \n',
		'0 : exit \n '
		'1 : Display number of peptides \n',
		'2 : Change lenght parameter \n',
		'3 : Change lag parameter \n',
		'4 : Distribution of lenghts\n',
		'5 : Distribution of aa frequency\n',
		'6 : Plot Mean Zscale\n',
		'7 : Plot ACC Zscale (ACP)\n',
		'8 : Tree of ACC euclidean distance\n',
		'9 : K-means of ACC euclidean distance \n')



	A = int(input ('->   '))


	#Display number of peptides
	if A == 1 :
		Nb_seq(BDs,Label)

	#Change lenght parameter
	if A == 2 :
		print( 'Actual lenght parameters : ',Inf,' to ',Sup)
		B = input ('\nChange ? (y/n) ')
		if B == 'y' :
			Inf = int(input('Minimal lenght : '))
			Sup = int(input('Maximal lenght : '))
			BDs,couleur,Label,Seq,Name = Init_DFs(Ps) # Reload initial sequences 
			BDs=All_Selec(Inf,Sup,BDs,Seq) 
			Crea_New(BDs,Label) 
			Run_ACC(BDs,Label,Seq,lag)
			Run_Z(Label,Seq)
			ACC = All_ACC(BDs,Label)
			Z =  All_Z(Label)
			MaJ_Bornes_Param(Inf,Sup)

	#Changer le lag
	if A == 3:
		print ('Actual lag : ',lag)
		lag = int(input ('Lag = '))
		Run_ACC(BDs,Label,Seq,lag)
		ACC = All_ACC(BDs,Label)
		MaJ_lag_Param(lag)

	#violinplot of lenght 
	if A == 4 :
		Aff_group(Label) 
		g = eval(input('Serval groups (ex : [0,1,2]) or just one (ex : 0) : '))
		if type(g) == int:
			LEN = [len(BDs[g].at[x,Seq[g]]) for x in range(len(BDs[g]))]
			sns.distplot(LEN, bins=20, rug=True)
			plt.title('Distibution lenght '+Label[g])
			plt.savefig('Out_File/Out_Length/Distri_Length_'+Label[g]+'.svg')
			plt.close()
			print('The figure of lenght distribution has been recorded: Out_File/Out_Length/Distri_Lenght_'+Label[g]+'.svg')
		else:
			violinplot_len(g,BDs,Seq,Label,couleur)

	#Frequency of aa
	if A== 5 :
		q=int(input('Mean frequency (0) By position (1 ) :'))
		Aff_group(Label)
		if q == 1 :
			g = int(input('Groups (ex : 0) : '))
			m = float(input('Frequency max (ex : 0.5) :'))
			Dist_Heat(g,m,BDs,Seq,Label)

		else :		
			h=eval(input('HeatMap (0) or Histogram (1) :'))
			if h==0 :
				g = eval(input('Groups (ex : [0,1,2]) : '))
				m = float(input('Frequency max (ex : 0.5) :'))
				Dist_Heat_Moy(g,m,BDs,Seq,Label)
			else :
				g = eval(input('Groups (ex : [0,1,2]) : '))
				plot_dist_feq_aa(g,BDs,Label,Seq,couleur)

	#Z-scale    
	if A == 6 :
		c = int(input ('Number dimension for represent mean z-scale 1 or 2 (use ACP) or 3  : '))
		Aff_group(Label) 
		g = eval(input('Groups (ex : [0,1,2]) : '))
		if c == 1:
			plot_Zscale_1D(g,Label,couleur,Z)
		if c == 2:
			b=int(input('Basic PCA (0) Edit parameter (1) : '))
			if b == 1 :
				Indmod=eval(input('Other groups for learn model ? No : 0 else : [0,1,2] : '))
				L_G=eval(input('Groups with bigger markers ? (oui : [0,1,2] / non : 0)'))
				n=int(input('Display name of first group ? (0:no 1:yes)'))
				air=int(input('Convexe area ? (0:non 1:oui)'))
				inarea=eval(input('File indicate if a sequence is in or out a convexe area in PCA representation ? No (0) Yes ([index of test group, index of area group]) :'))
				PCA_2D(Indmod,g, Z,BDs,Label,couleur,Name,Seq,n,air,'all','all',L_G,0,inarea)
			else :
				PCA_2D(0,g, Z,BDs,Label,couleur,Name,Seq,0,1,'all','all',0,0)
		if c == 3:
			plot_Zscale(g,Label,couleur,Z)

	#ACC ACP
	if A == 7:
		print ('Actual lag: ',lag)
		c = int(input ('Number of dimension 2  ou 3 : '))
		Aff_group(Label) 
		g = eval(input('Groups (ex : [0,1,2]) : '))
		if c == 2:
			b=int(input('Basic PCA (0) Edit parameter (1) : '))
			if b == 1 :
				Indmod=eval(input('Other groups for learn model ? No : 0 else : [0,1,2] : '))
				L_G=eval(input('Groups with bigger markers ? (oui : [0,1,2] / non : 0)'))
				z=input('Wich z-scale factor ? all : z1, z2 and z3  1 : z1  2 : z2  3 : z3  :')
				l=input('Which lag ? all : lag1, lag2, lag3 and lag4  1 : lag1  2 : lag2  3 : lag3  4 : lag4 :')
				auto=int(input('Only auto-covariance ? (yes : 1)'))
				n=int(input('Display name of first group ? (0:no 1:yes)'))
				air=int(input('Convexe area ? (0:non 1:oui)'))
				inarea=eval(input('File indicate if a sequence is in or out a convexe area in PCA representation ? No (0) Yes ([index of test group, index of area group]) :'))
				PCA_2D(Indmod,g, ACC,BDs,Label,couleur,Name,Seq,n,air,z,l,L_G,auto,inarea)
			else :
				PCA_2D(0,g, ACC,BDs,Label,couleur,Name,Seq,0,1,'all','all',0,0)

		if c == 3:
			PCA_3D(g, ACC,BDs,Label,couleur,Name,Seq)

	#Newick Tree
	if A == 8:
		Aff_group(Label) 
		g = eval(input('Group (ex : 0) : '))
		f = eval(input('Choose feature ACC or Mean of Z-scale ? (ACC or Z) : '))
		P = eval(input('Use 2D coordinate of a PCA ? (0 : no or 1 : yes)'))
		if P == 1 :
			Indmod=eval(input('Groups for learn model ? [0,1,2] : '))
			f=Use_PCA_Coordinate_as_Feature(Indmod,f,BDs)
		Make_newick_Tree(g,f,BDs,Label,Name)

	#K-Means
	if A == 9:
		Aff_group(Label) 
		g = eval(input('Groups (ex : [0,1,2]) : '))
		f = eval(input('Choose feature ACC or Mean of Z-scale ? (ACC or Z) : '))		
		P = eval(input('Use 2D coordinate of a PCA ? (0 : no or 1 : yes)'))
		if P == 1 :
			Indmod=eval(input('Groups for learn model ? [0,1,2] : '))
			f=Use_PCA_Coordinate_as_Feature(Indmod,f,BDs)
		k = eval(input('Number of cluster (ex : 4): '))
		Process_Kmean(g,f,BDs,Label,Name,couleur,k)



		

