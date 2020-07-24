# -*- coding: utf-8 -*-
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import researchpy as rp
import statsmodels.api as sm
from statsmodels.formula.api import ols

Label=['bSP','eSP','tSP','Group1','Cyclotide','Group2','mTP','cTP','Random']
#Label=['mTP','mTPHsapiens','mTPScerevisae','mTPAthaliana','mTPCreinhardtii','cTP','cTPAthaliana','cTPCreinhardtii','Random']
#Label=['mTPHsapiens','mTPScerevisae','mTPAthaliana','mTPCreinhardtii','cTPAthaliana','cTPCreinhardtii','Random']



#Titre=['bSP','eSP','tSP','Class II\nHA-RAMP','globular\nAMP','Class I\nHA-RAMP','mTP','cTP','Random']
#Titre=['bSP','eSP','tSP','Class IIHA-RAMP','globularAMP','Class IHA-RAMP','mTP','cTP','Random']
#Titre=['SP','SP','SP','Class II HA-RAMP','globular AMP','Class I HA-RAMP','TP','TP','Random']
Titre=['SP','SP','SP','Class II\nHA-RAMP','globular\nAMP','Class I\nHA-RAMP','mTP','cTP','Random']
#Titre=['SP','SP','SP','Class II\nHA-RAMP','globular\nAMP','Class I\nHA-RAMP','TP','TP','Random']
#Titre=['mTP','mTP\nH.sapiens','mTP\nS.cerevisae','mTP\nA.thaliana','mTP\nC.reinhardtii','cTP','cTP\nA.thaliana','cTP\nC.reinhardtii','Random']
#Titre=['mTP','mTPH.sapiens','mTPS.cerevisae','mTPA.thaliana','mTPC.reinhardtii','cTP','cTPA.thaliana','cTPC.reinhardtii','Random']
#Titre=['mTP\nH.sapiens','mTP\nS.cerevisae','mTP\nA.thaliana','mTP\nC.reinhardtii','cTP\nA.thaliana','cTP\nC.reinhardtii','Random']
#Titre=['mTPH.sapiens','mTPS.cerevisae','mTPA.thaliana','mTPC.reinhardtii','cTPA.thaliana','cTPC.reinhardtii','Random']

Seq=['Sequence','Sequence','Sequence','Sequence','Sequence','Sequence','Peptide','Peptide','Sequence']
#Seq=['Peptide','Peptide','Peptide','Peptide','Peptide','Peptide','Peptide','Peptide','Sequence']
#Seq=['Peptide','Peptide','Peptide','Peptide','Peptide','Peptide','Sequence']

couleur=['#390c6f','#0c6f54','#d6e51e','#7e1717','#e80ddb','#3274a1','#F79646','#9BBB59','#f1c10d']
# couleur=['#390c6f','#0c6f54','#d6e51e','#7e1717','#e80ddb','#3274a1','#F79646','#9BBB59','#f1c10d']
# couleur=['#0c6f54','#0c6f54','#0c6f54','#7e1717','#e80ddb','#3274a1','#F79646','#9BBB59','#f1c10d']
#couleur=['#a14000','#ff4312','#712424','#c90a12','#576a1c','#2c9838','#f1c10d']

Color={}
for i in range(len(couleur)):
	Color[Titre[i]]=couleur[i]

All_data=pd.DataFrame(columns=['Label','Proportion amphipatic helix','Hydrophobicity','Hydrique moment','Peptide length','net charge','number aa in amphipatic helix','helix net charge'])

for i in range(len(Label)):
	lab=Label[i]
	data=pd.read_csv('NewThreeHelix_'+lab+'.csv',sep=',')
	sour=pd.read_csv('../New_csv/'+lab+'.csv',sep='\t')

	for j in range (len(data)):
		seq=sour.at[j,Seq[i]]
		All_data.loc[len(All_data)]=[Titre[i],data.at[j,'PaaHheliquest'],data.at[j,'Hydro'],data.at[j,'HydrMom'],float(len(seq)),float(seq.count('R')+seq.count('K')-seq.count('D')-seq.count('E')),data.at[j,'PaaHheliquest']*float(len(seq)),data.at[j,'Hnetcharge']]

All_data.to_csv('All_data_for_figure2.csv',sep='\t')
# All_data.to_csv('All_data_for_figure3.csv',sep='\t')

