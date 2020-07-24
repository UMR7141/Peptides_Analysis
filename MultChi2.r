
Tab.Kmean <- read.csv("Out_File/Out_Kmean/Tab_Kmean.csv", sep = ",",header=TRUE,row.names=1)
#Tab.Kmean <- read.csv('/data/garrido/Figures/1ARTICLE/PrepFig/NewPrepFig020719/ValidFigure/Final Figure/New_Figure_ClassAB/FinalFig/Principal_program_Article/Out_File/Out_Kmean/Resultat/Article/Tab_Kmean.csv', sep = ',',header=TRUE,row.names=1)
Tab.Kmean <-(Tab.Kmean[c(1:dim(Tab.Kmean)[1]-1),c(1:dim(Tab.Kmean)[2]-1)])



print(chisq.test(Tab.Kmean))

# chi2.cal<-sum((Tab.Kmean-chisq.test(Tab.Kmean)$expected)^2/chisq.test(Tab.Kmean)$expected)
# p.val<-(round(suppressWarnings(chisq.test(Tab.Kmean,correct=TRUE)$p.value),Inf))
# print(p.val)
MultChi <- function(Tab.Kmean,Out_File){

k<-nrow(Tab.Kmean)
c<-ncol(Tab.Kmean)
class1.G<-NULL
class2.G<-NULL
class1.F<-NULL
test<-NULL
p.value<-NULL

for(col1 in 1:(c-1)){
	for(col2 in (col1+1):c){
		for (row1 in 1:(k)) {
			class1.G<-c(class1.G,dimnames(Tab.Kmean)[[2]][col1])
			class2.G<-c(class2.G,dimnames(Tab.Kmean)[[2]][col2])
			class1.F<-c(class1.F,dimnames(Tab.Kmean)[[1]][row1])
			sum1<-sum(Tab.Kmean[,col1])-Tab.Kmean[row1,col1]
			sum2<-sum(Tab.Kmean[,col2])-Tab.Kmean[row1,col2]
			ob<-matrix(c(Tab.Kmean[row1,col1],Tab.Kmean[row1,col2],sum1,sum2),2,2,byrow=TRUE)
			th<-suppressWarnings(chisq.test(ob)$expected)
			n.cases.inf.5<-length(which(th<5))
			if (n.cases.inf.5==0){test<-c(test,"Chi2")
				p.value<-c(p.value,round(suppressWarnings(chisq.test(ob,correct=TRUE)$p.value),Inf))}
			else{
				test<-c(test,"Fisher.exact")
				p.value<-c(p.value,round(fisher.test(ob)$p.value,9))}}}}
p.mult<-data.frame(class1.G,class2.G,class1.F,test,p.value)

p.mult2<-transform(p.mult,Holm= format(p.adjust(p.mult$p.value),scientific=TRUE,digits=3))


write.csv(p.mult2, file = Out_File,row.names=FALSE)
}

MultChi(Tab.Kmean,'Out_File/Out_Kmean/MultipleChi2byGroup.csv')
MultChi(t(Tab.Kmean),'Out_File/Out_Kmean/MultipleChi2byCluster.csv')