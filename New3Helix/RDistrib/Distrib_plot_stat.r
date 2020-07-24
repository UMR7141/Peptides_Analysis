library(ggplot2)

#The All_data files were generated with the python script ../Distrib.py
data <- read.csv('All_data_for_figure2.csv', sep = '\t')
# data <- read.csv('All_data_for_figure3.csv', sep = '\t')

couleur<-c('#0c6f54','#7e1717','#e80ddb','#3274a1','#F79646','#9BBB59','#f1c10d')
# couleur<-c('#a14000','#ff4312','#712424','#c90a12','#576a1c','#2c9838','#f1c10d')


data$Label <- factor(data$Label, levels = c("SP", "Class II\nHA-RAMP", "globular\nAMP","Class I\nHA-RAMP", "mTP", "cTP","Random"))
# data$Label <- factor(data$Label, levels = c('mTP\nH.sapiens','mTP\nS.cerevisae','mTP\nA.thaliana','mTP\nC.reinhardtii','cTP\nA.thaliana','cTP\nC.reinhardtii','Random'))

sink('statfile.txt')


a <- ggplot(data ,aes(x=Label, y=Peptide.length, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))
a + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)


print('Peptide.length')

dd  <- pairwise.wilcox.test(data$Peptide.length, 
                           data$Label,p.adjust.method ="holm" )
print(dd)




b <- ggplot(data ,aes(x=Label, y=Hydrophobicity, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))
b + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)


print('Hydrophobicity')

dd  <- pairwise.wilcox.test(data$Hydrophobicity, 
                           data$Label,p.adjust.method ="holm" )
print(dd)




c <- ggplot(data ,aes(x=Label, y=net.charge, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))
c + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)



print('net.charge')

dd  <- pairwise.wilcox.test(data$net.charge, 
                           data$Label,p.adjust.method ="holm" )

print(dd)



d <- ggplot(data ,aes(x=Label, y=Hydrophobique.moment, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))
d + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)


print('Hydrophobique.moment')

dd  <- pairwise.wilcox.test(data$Hydrophobique.moment, 
                           data$Label,p.adjust.method ="holm" )

print(dd)


print('SP')
print(nrow(data[data$Label == 'SP', ]))
print('Class II\nHA-RAMP')
print(nrow(data[data$Label == 'Class II\nHA-RAMP', ]))
print('globular\nAMP')
print(nrow(data[data$Label == 'globular\nAMP', ]))
print('Class I\nHA-RAMP')
print(nrow(data[data$Label == 'Class I\nHA-RAMP', ]))
print('mTP')
print(nrow(data[data$Label == 'mTP', ]))
print('cTP')
print(nrow(data[data$Label == 'cTP', ]))
print('Random')
print(nrow(data[data$Label == 'Random', ]))


# print('mTP  H.sapiens')
# print(nrow(data[data$Label == 'mTP\nH.sapiens', ]))
# print('mTP  S.cerevisae')
# print(nrow(data[data$Label == 'mTP\nS.cerevisae', ]))
# print('mTP  A.thaliana')
# print(nrow(data[data$Label == 'mTP\nA.thaliana', ]))
# print('mTP  C.reinhardtii')
# print(nrow(data[data$Label == 'mTP\nC.reinhardtii', ]))
# print('cTP  A.thaliana')
# print(nrow(data[data$Label == 'cTP\nA.thaliana', ]))
# print('cTP  C.reinhardtii')
# print(nrow(data[data$Label == 'cTP\nC.reinhardtii', ]))
# print('Random')
# print(nrow(data[data$Label == 'Random', ]))

data <- data[data$number.aa.in.amphipatic.helix != 0, ]

options(repr.plot.width = 40, repr.plot.height = 8)

bp <- ggplot(data ,aes(x=Label, y=number.aa.in.amphipatic.helix, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))
bp + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)



print('number.aa.in.amphipatic.helix')

dd  <- pairwise.wilcox.test(data$number.aa.in.amphipatic.helix, 
                           data$Label,p.adjust.method ="holm" )

print(dd)






bp <- ggplot(data ,aes(x=Label, y=Proportion.amphipatic.helix, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))
bp + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)

print('Proportion.amphipatic.helix')

dd  <- pairwise.wilcox.test(data$Proportion.amphipatic.helix, 
                           data$Label,p.adjust.method ="holm" )

print(dd)




bp <- ggplot(data ,aes(x=Label, y=helix.net.charge, colour=Label,fill=Label))+
    geom_jitter(width=0.25)+
    geom_boxplot(alpha=0.25, outlier.alpha=0, notch=TRUE) +
    stat_summary(fun.y=mean, colour="black", geom="point", 
                 shape=18, size=3,show.legend  = FALSE)+
    theme_classic()+
    theme(legend.position="none")+
    theme(axis.text = element_text(angle=0,size=12),
        axis.title=element_text(size=14,face="bold"))

bp + scale_fill_manual(values=couleur) + scale_color_manual(values=couleur)


print('helix.net.charge')

dd  <- pairwise.wilcox.test(data$helix.net.charge, 
                           data$Label,p.adjust.method ="holm" )

print(dd)


# ggsave(file="test.svg", plot=image, width=10, height=8)




print('SP')
print(nrow(data[data$Label == 'SP', ]))
print('Class II\nHA-RAMP')
print(nrow(data[data$Label == 'Class II\nHA-RAMP', ]))
print('globular\nAMP')
print(nrow(data[data$Label == 'globular\nAMP', ]))
print('Class I\nHA-RAMP')
print(nrow(data[data$Label == 'Class I\nHA-RAMP', ]))
print('mTP')
print(nrow(data[data$Label == 'mTP', ]))
print('cTP')
print(nrow(data[data$Label == 'cTP', ]))
print('Random')
print(nrow(data[data$Label == 'Random', ]))

# print('mTP  H.sapiens')
# print(nrow(data[data$Label == 'mTP\nH.sapiens', ]))
# print('mTP  S.cerevisae')
# print(nrow(data[data$Label == 'mTP\nS.cerevisae', ]))
# print('mTP  A.thaliana')
# print(nrow(data[data$Label == 'mTP\nA.thaliana', ]))
# print('mTP  C.reinhardtii')
# print(nrow(data[data$Label == 'mTP\nC.reinhardtii', ]))
# print('cTP  A.thaliana')
# print(nrow(data[data$Label == 'cTP\nA.thaliana', ]))
# print('cTP  C.reinhardtii')
# print(nrow(data[data$Label == 'cTP\nC.reinhardtii', ]))
# print('Random')
# print(nrow(data[data$Label == 'Random', ]))



















