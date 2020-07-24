#Command line : Rscript --vanilla acc.R -f <csv-file> -a <column-name> -l <lag> -o <out-file>
#!/usr/bin/env Rscript
#install.packages("optparse")
#install.package('protr')
#install.packages("stringr", dependencies=TRUE)
#install.packages("scatterplot3d")
library('protr')
library('stringr')
library("scatterplot3d")
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="csv file", metavar="character"),
  make_option(c("-a", "--accession"), type="character", default=NULL, 
              help='name column sequence', metavar="character"),
  make_option(c("-l", "--lag"), type="integer", default=NULL,
              help="l refers to lag, which is the interval between residues being compared"),
  make_option(c("-o", "--out-file"), type="character", default=NULL,
              help="out csv file")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


data <- read.csv(opt$f, sep = '\t')

mat_vect=c()

for (s in (1:dim(data)[1]))
	{
	seq = as.character(data[s,opt$a])
	mat = rbind(AAindex[390,str_sub(seq,1,1)],AAindex[391,str_sub(seq,1,1)],AAindex[392,str_sub(seq,1,1)])
	for (i in (2:nchar(seq)))
		{	
		if (str_sub(seq,i,i)=="Z"|str_sub(seq,i,i)=="U"|str_sub(seq,i,i)=="O"|str_sub(seq,i,i)=="J"|str_sub(seq,i,i)=="B"|str_sub(seq,i,i)=="X" )
			{
			cat('Warning : unrecognized amino acid : ',str_sub(seq,i,i),'\n')
			}
		
		mat = cbind(mat,rbind(AAindex[390,str_sub(seq,i,i)],AAindex[391,str_sub(seq,i,i)],AAindex[392,str_sub(seq,i,i)]))
		}
	mat = t(as.matrix(acc(mat,opt$l)))
	if (length(mat_vect)==0)
		{
		mat_vect = mat
		}
	else
		{
		mat_vect = rbind(mat_vect,mat)
		}
	}

# Change name of columns
for (x in (1:length(colnames(mat_vect)))){colnames(mat_vect)[x]<-sub('scl','z',colnames(mat_vect)[x])} 
write.csv(mat_vect, file = opt$o,row.names=FALSE)

