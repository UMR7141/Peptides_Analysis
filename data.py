# -*- coding: utf-8 -*-

# You can add your data just add P_** as is describe lower and add this in Ps list 

# P_** = [Path to csv table, Separator, Desired color, Desired label, Name of column contained sequences,  Name of column contained names of sequences]
# P_** = [Path, Separator, Color, Label, Index_Seq, Index_Name]

# Ps = [P_1,P_2,...,P_n] #List of used data for analyze


################
# Initial Data #
################



#HA-RAMP Helical-Amphiphilic Ribosomal-Associated Antimicrobial Peptide

P_HARAMPClassI = ['In_File/HA-RAMP/Group2.csv','\t','#3274a1','HA-RAMPclassI','Sequence','Index']
P_HARAMPClassII = ['In_File/HA-RAMP/Group1.csv','\t','#ec513f','HA-RAMPclassII','Sequence','Index']
P_Mag=['In_File/HA-RAMP/Magainin.csv','\t','#6372dc','Magainin','Sequence','Index']
P_Cer=['In_File/HA-RAMP/Cecropin.csv','\t','#ea722a','Cecropin','Sequence','Index']
P_Der=['In_File/HA-RAMP/Dermaseptin.csv','\t','#efe80f','Dermaseptin','Sequence','Index']
P_Esc1=['In_File/HA-RAMP/Esculentin-1.csv','\t','#bdc2c2','Esculentin-1','Sequence','Index']
P_Esc2=['In_File/HA-RAMP/Esculentin-2.csv','\t','#bdc2c2','Esculentin-2','Sequence','Index']
P_Max=['In_File/HA-RAMP/Maximin.csv','\t','#d61625','Maximin','Sequence','Index']
P_Def=['In_File/HA-RAMP/Defensin.csv','\t','#0f23f0','Defensin','Sequence','Index']
P_Bac=['In_File/HA-RAMP/Bacteriocin IIa.csv','\t','#7ee0dd','BacteriocinIIa','Sequence','Index']
P_Lat=['In_File/HA-RAMP/Latracin.csv','\t','#43a741','Latracin','Sequence','Index']
P_Mas=['In_File/HA-RAMP/Mastoparan.csv','\t','#7ee081','Mastoparan','Sequence','Index']
P_Mori=['In_File/HA-RAMP/Moricin.csv','\t','#4e5252','Moricin','Sequence','Index']
P_Cat=['In_File/HA-RAMP/Cathelicidin.csv','\t','#e05cd8','Cathelicidin','Sequence','Index']
P_Oce=['In_File/HA-RAMP/Ocellatin.csv','\t','#ee2168','Ocellatin','Sequence','Index']
P_Pleu=['In_File/HA-RAMP/Pleurocidin.csv','\t','#67e3e7','Pleurocidin','Sequence','Index']
P_Rana=['In_File/HA-RAMP/Ranatuerin.csv','\t','#9eec3f','Ranaturien','Sequence','Index']
P_Brev2=['In_File/HA-RAMP/Brevinin-2.csv','\t','#e07ed7','Brevinin-2','Sequence','Index']
P_Aur=['In_File/HA-RAMP/Aurein.csv','\t','#eb8218','Aurein','Sequence','Index']
P_Cae=['In_File/HA-RAMP/Caerin.csv','\t','#ec863c','Caerin','Sequence','Index']
P_Gae=['In_File/HA-RAMP/Gaegurin.csv','\t','#9be7d3','Gaegurin','Sequence','Index']
P_His=['In_File/HA-RAMP/Histatin.csv','\t','#ec1f15','Histatin','Sequence','Index']
P_Phy=['In_File/HA-RAMP/Phylloseptin.csv','\t','#efe80f','Phyloseptin','Sequence','Index']
P_Brev1=['In_File/HA-RAMP/Brevinin-1.csv','\t','#26f389','Brevinin-1','Sequence','Index']
P_Bom=['In_File/HA-RAMP/Bombinin.csv','\t','#dccc63','Bombinin','Sequence','Index']
P_Mel=['In_File/HA-RAMP/Melittin.csv','\t','#f1150e','Melittin','Sequence','Index']
P_Cla=['In_File/HA-RAMP/Clavanin.csv','\t','#9fc17c','Clavanin','Sequence','Index']
P_Dermi=['In_File/HA-RAMP/Dermicidin.csv','\t','#ed07f5','Dermicidin','Sequence','Index']
P_Pse=['In_File/HA-RAMP/Pseudin.csv','\t','#956797','Pseudin','Sequence','Index']
P_Thi=['In_File/HA-RAMP/Thionin.csv','\t','#679497','Thionin','Sequence','Index']
P_Upe=['In_File/HA-RAMP/Uperin.csv','\t','#679497','Uperin','Sequence','Index']
P_Par=['In_File/HA-RAMP/Pardaxin.csv','\t','#679497','Pardaxin','Sequence','Index']
P_Tra=['In_File/HA-RAMP/Transferrin.csv','\t','#679497','Transferrin','Sequence','Index']

#TP Targeting Peptide

P_TP=['In_File/TP/TP.csv','\t','#000000','TP','Peptide','Index']
P_cTP=['In_File/TP/cTP.csv','\t','#6dc068','cTP','Peptide','Index']
P_cAra=['In_File/TP/cTPArabidopsis.csv','\t','#6bcd68','cTPAthaliana','Peptide','Index']
P_cChla=['In_File/TP/cTPChlamydomonas.csv','\t','#1d9219','cTPCreinhardtii','Peptide','Index']
P_mTP=['In_File/TP/mTP.csv','\t','#ff902e','mTP','Peptide','Index']
P_mChla=['In_File/TP/mTPChlamydomonas.csv','\t','#ad1616','mTPCreinhardtii','Peptide','Index']
P_mAra=['In_File/TP/mTPArabidopsis.csv','\t','#d46262','mTPAthaliana','Peptide','Index']
P_mHom=['In_File/TP/mTPHomo.csv','\t','#e55018','mTPHsapiens','Peptide','Index']
P_mCer=['In_File/TP/mTPCerevisae.csv','\t','#f89136','mTPScerevisae','Peptide','Index']


#SP Signal Peptide

P_SP=['In_File/SP/SP.csv','\t','#000000','SP','Sequence','Index']
P_bSP=['In_File/SP/bSP.csv','\t','#390c6f','bSP','Sequence','Index']
P_tSP=['In_File/SP/tSP.csv','\t','#d6e51e','tSP','Sequence', 'Index']
P_eSP=['In_File/SP/eSP.csv','\t','#0c6f54','eSP','Sequence','Index']


#Other

P_gAMP=['In_File/Other/Cyclotide.csv','\t','#e80ddb','Cyclotide','Sequence','Index']
P_Random=['In_File/Other/Random.csv',',','#d6a4d6','Random','Sequence','Index']


#List of data using in programm

Ps=[P_HARAMPClassI,P_HARAMPClassII,P_Mag,P_Cer,P_Der,P_Esc1,P_Esc2,P_Max,P_Def,P_Bac,P_Lat,P_Mas,P_Def,P_Bac,P_Lat,P_Mas,P_Mori,P_Cat,P_Oce,P_Pleu,P_Rana,P_Brev2,P_Aur,P_Cae,P_Gae,P_His,P_Phy,P_Brev1,P_Bom,P_Mel,P_Cla,P_Dermi,P_Pse,P_Thi,P_Upe,P_Par,P_Tra,P_TP,P_cTP,P_cAra,P_cChla,P_mTP,P_mAra,P_mHom,P_mCer,P_SP,P_bSP,P_tSP,P_eSP,P_gAMP,P_Random]