#Command line : Rscript --vanilla Z_scale_Wold.R -f <csv-file> -a <column-name> -o <out-file>
#!/usr/bin/env Rscript
#install.packages("optparse")
#install.package('protr')
#install.packages("stringr", dependencies=TRUE)
library('protr')
library('stringr')
library("optparse")



option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="csv file", metavar="character"),
  make_option(c("-a", "--accession"), type="character", default=NULL, 
              help='name column sequence', metavar="character"),
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
		mat = cbind(mat,rbind(AAindex[390,str_sub(seq,i,i)],AAindex[391,str_sub(seq,i,i)],AAindex[392,str_sub(seq,i,i)]))
		}
	mat = rowMeans(mat)

	if (length(mat_vect) == 0)
		{
		mat_vect = mat
		}
	else
		{
		mat_vect = rbind(mat_vect,mat)
		}
	}

write.csv(mat_vect, file = opt$o,row.names=FALSE)

