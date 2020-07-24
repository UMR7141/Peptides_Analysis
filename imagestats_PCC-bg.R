#calculate Pearson Correllation Coefficients for Venus compared to chlorophyll and mitotracker channels in batch
library(tiff)
library(imager)
source("/Users/Oliver/Google Drive/Postdoc/Paris/R scripts/barchart.oli.R")
source("/Users/Oliver/R/JTS-func.R")
source("/Users/Oliver/R/SE.R")

#scriptsettings
substractVenus <- TRUE
substractChlorophyll <- TRUE
substractMitoTracker <- TRUE
SDmult <- 1

#substractionfile
backgroundsvalues <- read.csv("/Users/Oliver/Google Drive/Postdoc/Paris/Publications/1 - Bioinfo + M2/cells/BackgroundSubstractions.csv",header = TRUE)

#filepaths
path_611 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/pMO611 - Venus/stats/"
#path_611man <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/pMO611 - Venus/stats -bg man nocp/"
path_017 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/pODC17 - AcTP/stats/"
path_034 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/ AcTP+23 (pODC34)/stats/"
path_076 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/cTP cleavage site (pODC76)/stats/"
path_055 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/pODC55 - CAG2+23/stats/"
#path_009 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/controls/pODC9 - no Venus/stats/"

path_101 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/AMP-cTP/M2-cTP (pODC101)/stats/"
path_121 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/AMP-cTP/B15-cTP (pODC121)/stats/"

path_186 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/AMP-cTP/B2I-cTP (pODC186)/stats images/"
path_187 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/AMP-cTP/S1D-cTP (pODC187)/stats images/"
path_188 <- "/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/AMP-cTP/EHF-cTP (pODC188)/stats images/"



paths <- apropos("path_") 

imagearrayslices <- c()
#loop through all paths
for (path.i in 1:length(paths)) {
  setwd(get(paths[path.i]))
  current.files <- list.files(pattern=".tif")
  
  current.images <- lapply(current.files, readTIFF, all=TRUE)
  #At this point, we have a list of all images in the directory (aka of this strain).
  #Each image is present as a list of length three, corresponding to the three channels.
  #Channel 1 is Venus 
  #Channel 2 is Chlorophyll
  #Channel 3 is Mitotracker.
  #Each channel contains the image data as a pixel matrix with pixel values between 0 and 65536.

  
  PCC_VM <- c()
  PCC_VC <- c()
  PCC_MC <- c()
  maxMTpx <- c()
  for (image.i in 1:length(current.images)) {
    test.image <- current.images[[image.i]]
    
    #clean up mysterious multiple arrays per image
    #readTIFF sometimes generates two image array slices and one slice that is all 0s
    #so I'm finding the one that has values >0 and assume this is the real image
    #bitsize <- rep(16,3)
    for (test.i in 1:3) {
      if (!is.na(dim(test.image[[test.i]])[3])) {
        if (any(test.image[[test.i]][,,1]>0)) {
          test.image[[test.i]] <- test.image[[test.i]][,,1]
        } else if (any(test.image[[test.i]][,,2]>0)) {
          test.image[[test.i]] <- test.image[[test.i]][,,2]
        } else if (any(test.image[[test.i]][,,3]>0)) {
          test.image[[test.i]] <- test.image[[test.i]][,,3]
        } else {
          imagearrayslices <- c(imagearrayslices, paste(path.i, image.i, test.i, sep = "_"))
        }
        
      }
     }
    

    #substract background
    if (any(c(substractVenus, substractChlorophyll, substractMitoTracker))) {
      
      #find correct line in backgroundvalues dataframe
      strain <- strsplit(paths[path.i],split = "_")[[1]][2]
      if (strain != "611man") {
        
      if (strsplit(strain, split = "")[[1]][1]=="0") strain <- paste(strsplit(strain, split = "")[[1]][-1], collapse = "")
      date <- paste0(strsplit(current.files[image.i],split = "")[[1]][1:8], collapse ="")
      date <- strsplit(date, split="-")[[1]]
      imageID <-  strsplit(current.files[image.i],split = "-")[[1]]
      imageID <- strsplit(imageID[length(imageID)], split=".tif")[[1]]
    
      thisBGline <- which(backgroundsvalues$strain==strain)
      thisBGline <- thisBGline[which(backgroundsvalues$date[thisBGline]==date)]
      thisBGline <- thisBGline[which(backgroundsvalues$image[thisBGline]==imageID)]
      
      if (length(thisBGline) != 1) {
        print(paste("BGline error at", strain, date, imageID))
        break
      }
      

      if (backgroundsvalues$bit[thisBGline]==8)
      { bitsizer <- 1/256 } else {bitsizer <- 1/65536}    
      #substract background values
      if (substractVenus) test.image[[1]] <- test.image[[1]] - backgroundsvalues$Venus.mean[thisBGline] * bitsizer - backgroundsvalues$Venus.sd[thisBGline]* bitsizer *SDmult
      if (substractChlorophyll) test.image[[2]] <- test.image[[2]] - backgroundsvalues$Chloroplast.mean[thisBGline]  * bitsizer    - backgroundsvalues$Chloroplast.sd[thisBGline]  * bitsizer *SDmult
      if (substractMitoTracker) test.image[[3]] <- test.image[[3]] - backgroundsvalues$MitoTracker.mean[thisBGline]* bitsizer - backgroundsvalues$MitoTracker.sd[thisBGline]* bitsizer*SDmult


      #reset negative values to zero
      for (channel.i in c(1:3)) {
        #test.image[[channel.i]] <- round(test.image[[channel.i]], digits=0)
        if (any(test.image[[channel.i]] <0 )) test.image[[channel.i]][which(test.image[[channel.i]]<0)] <- 0
      }
    }}
    maxMTpx[image.i] <- max(test.image[[3]])
    (PCC_VM[image.i] <- cor(x=as.vector(test.image[[3]]), y=as.vector(test.image[[1]]), method = "pearson"))
    (PCC_VC[image.i] <- cor(x=as.vector(test.image[[2]]), y=as.vector(test.image[[1]]), method = "pearson"))
    (PCC_MC[image.i] <- cor(x=as.vector(test.image[[2]]), y=as.vector(test.image[[3]]), method = "pearson"))
    if (any(c(is.na(PCC_VM[image.i]), is.na(PCC_VC[image.i]), is.na(PCC_MC[image.i])))) {
      print("error, PCC is NA")
      break
    }
  }
  
  names(maxMTpx) <- current.files
  names(PCC_VM) <- current.files
  names(PCC_VC) <- current.files
  names(PCC_MC) <- current.files
  names(current.images) <- current.files
  #store values in a variable named after the respective strain
  assign(x= paste0("maxMTpx_", strsplit(paths[path.i], split="_")[[1]][2]),value=maxMTpx)
  assign(x= paste0("PCC_VM_", strsplit(paths[path.i], split="_")[[1]][2]),value=PCC_VM)
  assign(x= paste0("PCC_VC_", strsplit(paths[path.i], split="_")[[1]][2]),value=PCC_VC)
  assign(x= paste0("PCC_MC_", strsplit(paths[path.i], split="_")[[1]][2]),value=PCC_MC)
  assign(x= paste0("images_", strsplit(paths[path.i], split="_")[[1]][2]), value = current.images)
}




order <- c("611", "017", "034", "076", "055", "121", "101", "186", "187", "188")
#order <- c("611man", "611")
labels <- c("Venus only", "RBCA-cTP", "RBCA-cTP+", "RBCA-cs", "CAG2-mTP+", "Bacillocin 1580", "Magainin II", "Brevinin-2ISb", "Sarcotoxin-1D", "Enterocin HF")
#labels <- order
means_VM <- c()
means_VC <- c()
means_MC <- c()
means <- c()
SEs_VM <- c()
SEs_VC <- c()
SEs_MC <- c()
SEs <- c()
long.labels <- c()
p.values <- c()
PCCs_VM <- c()
PCCs_VC <- c()
PCCs_MC <- c()
PCCs_names <- c()
for (i in 1:length(order)) {
  PCCs_names <- c(PCCs_names, rep(order[i], length(get(paste0("PCC_VM_", order[i])))))
  PCCs_VM <- c(PCCs_VM, get(paste0("PCC_VM_", order[i])))
  PCCs_VC <- c(PCCs_VC, get(paste0("PCC_VC_", order[i])))
  PCCs_MC <- c(PCCs_MC, get(paste0("PCC_MC_", order[i])))
  means_VM[i] <- mean(get(paste0("PCC_VM_", order[i])) )
  means_VC[i] <- mean(get(paste0("PCC_VC_", order[i])))
  means_MC[i] <- mean(get(paste0("PCC_MC_", order[i])))
  means <- c(means, means_VM[i], means_VC[i])
  SEs_VM[i] <- SE(get(paste0("PCC_VM_", order[i])))
  SEs_VC[i] <- SE(get(paste0("PCC_VC_", order[i])))
  SEs_MC[i] <- SE(get(paste0("PCC_MC_", order[i])))
  SEs <- c(SEs, SEs_VM[i], SEs_VC[i])
  long.labels <- c(long.labels, paste(order[i], "mT", sep=" -"), paste(labels[i], "chl", sep=" -"))
  test <- t.test(x=get(paste0("PCC_VM_", order[i])), y=get(paste0("PCC_VC_", order[i])), alternative="two.sided")
  p.values[i] <- test$p.value
}
spaces <- rep(c(1,0), length(labels))
colours <- rep(c("cyan", "magenta"), length(labels))

setwd("/Users/Oliver/Google Drive/Postdoc/Paris/data/Venus microscopy/presentables/stats images/")



PCCs <- data.frame(PCCs_names, PCCs_VM, PCCs_VC, PCCs_MC, cp=((PCCs_names=="034")+(PCCs_names=="121")), mt=((PCCs_names=="055")+(PCCs_names=="101")))


#ANOVA
anova_VM <- aov(PCCs_VM ~ PCCs_names, data = PCCs) 
anova_VC <- aov(PCCs_VC ~ PCCs_names, data = PCCs) 
anova_MC <- aov(PCCs_MC ~ PCCs_names, data = PCCs) 

Tukey_VM <- TukeyHSD(anova_VM, conf.level = 0.95)
Tukey_VM$PCCs_names[Tukey_VM$PCCs_names[,4]<0.05,4]*100
Tukey_VC <- TukeyHSD(anova_VC, conf.level = 0.95)
Tukey_VC$PCCs_names[Tukey_VC$PCCs_names[,4]<0.05,4]*100

ttest_VM <- t.test(PCCs_VM ~mt, data = PCCs)
ttest_VC <- t.test(PCCs_VC ~cp, data = PCCs)


#plot the data

pdf("20191127_PCCs-bgVCM_bychannel_butterfly_v2.pdf")

sigdig <- 2
lineoffset <- 0.25
xoffset <- 0.03

x_control <- 0.74
x_sample <- 0.84

ymax <- length(order)*1.2+0.5

y_strain <- ymax-c(1:length(order))*1.2+0.2

par(mfrow=c(1,1))
layout(matrix(c(1,2,3), ncol=3), widths = c(3,1.55,3))
par(mar=c(1,1,1,1),
    oma=c(6,1,1,1))

plotorder <- c("VM", "VC")
cex <- 1.5
for (plotnumber in 1:2) {
  thisplot <- plotorder[plotnumber]
  
  barchart.oli(means=(get(paste0("means_", thisplot)))[length(order):1], SE = get(paste0("SEs_", thisplot))[length(order):1],  labels="", ymin=0, ymax =ymax,xlim=c(c(1.1,0)[plotnumber],c(0,1.1)[plotnumber]), col=c("cyan", "magenta", "white")[plotnumber], las=3, horiz=TRUE, axes=FALSE)
  axis(1, cex=cex)
  points(y=c(-10,ymax-0.4), x=rep(mean(PCCs$PCCs_MC)+1.96*sd(PCCs$PCCs_MC),2), type="l", lty=2)
  #points(x=c(-10,100), y=rep(mean(PCCs$PCCs_MC),2), type="l", lty=1)
  points(y=c(-10,ymax-0.4), x=rep(mean(PCCs$PCCs_MC)-1.96*sd(PCCs$PCCs_MC),2), type="l", lty=2)
  
  text(y=ymax-0.2, x=0.4, labels=c("Venus vs. MitoTracker", "Venus vs. Chlorophyll")[plotnumber], cex=cex)
  
  
  #add p-values
  #controls <- c("055", "034")
  samples <- c("055_1", "034_2","101_1", "121_2", "186_1", "187_2", "188_2")
  

  #control.y <- y_strain[which(order==controls[plotnumber])]
  sample.y <- y_strain[which(order==samples[plotnumber])]
  names(y_strain) <- order
  
  #points(y=c(min(y_strain), max(y_strain)), x=rep(x_control,2), type="l", lty=1)
  #points(y=rep(y_strain[which(order == controls[plotnumber])],2), x=c(x_control,(get(paste0("means_", thisplot)))[which(order == controls[plotnumber])]), type="l", lty=1)
  
  addme <- 0
  for (sample.i in 1:length(samples)) {
    
    #plot vertical lines emanating from sample
    this.id <- strsplit(samples[sample.i], split = "_")[[1]] 
    if (as.numeric(this.id[2]) != plotnumber) next
    #points(y=c(min(y_strain), max(y_strain)), x=rep(x_sample,2), type="l", lty=1)
    points(y=rep(y_strain[which(order == this.id[1])],2), x=c(x_control+addme,(get(paste0("means_", thisplot)))[which(order == this.id[1])]), type="l", lty=1)
  
    
  #p.controls <- c()
  p.samples <- c()
  
  for (strain.i in 1:length(order)) {
    this.strain <- order[strain.i]
    
    #if (this.strain != controls[plotnumber]) {
    #  this.tukey <- get(paste0("Tukey_", plotorder[plotnumber]))
    #  control_lines <- grep(controls[plotnumber], rownames(this.tukey$PCCs_names))
    #  this.line <- control_lines[grep(this.strain, rownames(this.tukey$PCCs_names[control_lines,]))]
    #  this.p <- this.tukey$PCCs_names[this.line,"p adj"]
    #  p.controls[strain.i] <- this.p
    #  if (this.p < 0.05) {
    #    print.p <- formatC(signif(this.p, digits = sigdig))
    #    print.p <- strsplit(print.p, split = "e")[[1]]
    #    print.p <- bquote(.(print.p[1])%.%10^.(as.numeric(print.p[2])))
    #    
    #    #plot lines
    #    this.y <- y_strain[strain.i]
    #    if ((this.y>sample.y&&control.y<sample.y)) {
    #      points(y=c(this.y, sample.y+lineoffset), x=rep(x_control,2), type="l", lty=1)
    #      points(y=c(sample.y-lineoffset, control.y), x=rep(x_control,2), type="l", lty=1)
    #      points(y=c(sample.y+lineoffset,sample.y-lineoffset), x=rep(x_control,2), type="l", lty=3)
    #    } else if (this.y<sample.y&&control.y>sample.y) {
    #      points(y=c(this.y, sample.y-lineoffset), x=rep(x_control,2), type="l", lty=1)
    #      points(y=c(sample.y+lineoffset, control.y), x=rep(x_control,2), type="l", lty=1)
    #      points(y=c(sample.y+lineoffset,sample.y-lineoffset), x=rep(x_control,2), type="l", lty=3)
    #    } else {
    #      points(y=c(this.y, control.y), x=rep(x_control,2), type="l", lty=1)
    #    }
    #    points(y=rep(this.y,2), x=c(x_control,x_control-0.02), type="l", lty=1)
    #    if (this.p > 0) {
    #      #text(x= x_control-0.05, y= y_strain[strain.i],labels = print.p, srt=c(90, 270)[plotnumber])
    #    } else {
    #      #text(x= x_control-0.05, y= y_strain[strain.i],labels = bquote("<"*.(signif(2.62, digits=sigdig))%.%10^-14), srt=c(90, 270)[plotnumber])
    #    }
    #  }
    #}
    
    if (this.strain != this.id[1]) {
      this.tukey <- get(paste0("Tukey_", plotorder[plotnumber]))
      sample_lines <- grep(this.id[1], rownames(this.tukey$PCCs_names))
      this.line <- sample_lines[grep(this.strain, rownames(this.tukey$PCCs_names[sample_lines,]))]
      this.p <- this.tukey$PCCs_names[this.line,"p adj"]
      p.samples[strain.i] <- this.p
      if (this.p < 0.05) {
        #print.p <- formatC(signif(this.p, digits = sigdig))
        #print.p <- strsplit(print.p, split = "e")[[1]]
        #print.p <- bquote(.(print.p[1])%.%10^.(as.numeric(print.p[2])))
        #points(y=c(y_strain[which(order==this.id[1])], y_strain[strain.i]), x=rep(x_sample,2), type="l", lty=1)
        #points(y=rep(y_strain[strain.i],2), x=c(x_sample,x_sample-0.02), type="l", lty=1)
        #if (this.p > 0) {
          #text(x= x_sample-0.05, y= y_strain[strain.i],labels = print.p, srt=c(90, 270)[plotnumber])
        #} else {
          #text(x= x_sample-0.05, y= y_strain[strain.i],labels = bquote("<"*.(signif(2.62, digits=sigdig))%.%10^-14), srt=c(90, 270)[plotnumber])
        #}
        
            #plot lines
            current.samples <- samples[grep(paste0("_", plotnumber),samples)]
            crossed.samples <- current.samples[-(1:grep(this.id[1], current.samples))]
            if (length(grep(this.strain, crossed.samples))>0) crossed.samples <- crossed.samples[-grep(this.strain, crossed.samples)]
            
            
            #check if any lines are crossed
            if (length(crossed.samples)>0) {
              crossme.y <- c()
              for (crossers in 1:length(crossed.samples)) {
                this.crossed.sample <- strsplit(crossed.samples[crossers], split="_")[[1]][1]
                crossme.y[crossers] <- y_strain[which(order==this.crossed.sample)]
              }
        
              this.y <- y_strain[strain.i]
              
              passed.ys <- y_strain[strain.i:which(order==this.id[1])]
              
              crossed.ys <- c()
              for (crosses in 1:length(crossme.y)) {
                if (length(which(passed.ys == crossme.y[crosses]))>0) crossed.ys <- c(crossed.ys, crossme.y[crosses])
              }
            } else {crossed.ys <- c()}
              
            if (length(crossed.ys) > 0) { 
              #plot top line
              top.y <- max(y_strain[c(strain.i, which(order==this.id[1]))])
              points(y=c(top.y, max(crossed.ys)+lineoffset), x=rep(x_control+addme,2), type="l", lty=1)
              
              #plot bottom line
              bot.y <- min(y_strain[c(strain.i, which(order==this.id[1]))])
              points(y=c(bot.y, min(crossed.ys)-lineoffset), x=rep(x_control+addme,2), type="l", lty=1)
              
              #plot intermediate lines
              if (length(crossed.ys)>1) {
              for (linecounter in 1:(length(crossed.ys)-1)) {
                top.y <- max(crossed.ys[c(linecounter, linecounter+1)])
                bot.y <- min(crossed.ys[c(linecounter, linecounter+1)])
                points(y=c(bot.y+lineoffset, top.y-lineoffset), x=rep(x_control+addme,2), type="l", lty=1)
              }}
              
              #plot connecting lines
              for (linecounter in 1:length(crossed.ys)) {
                points(y=c(crossed.ys[linecounter]+lineoffset, crossed.ys[linecounter]-lineoffset), x=rep(x_control+addme,2), type="l", lty=3)
              }
              
              
            
            } else {
              points(y=c(y_strain[which(order==this.id[1])], y_strain[strain.i]), x=rep(x_control+addme,2), type="l", lty=1)
            }
            #plot indicators
            points(y=rep(y_strain[strain.i],2), x=c(x_control+addme,x_control+addme-0.02), type="l", lty=1)
      }
    }
  }
  
  #print.p <- formatC(signif(max(p.controls[which(p.controls<0.05)], na.rm = TRUE), digits = sigdig))
  #if (as.numeric(print.p)>0) {
  #  print.p <- strsplit(print.p, split = "e")[[1]]
  #  print.p <- bquote(italic(p)-values<=.(print.p[1])%.%10^.(as.numeric(print.p[2])))
  #} else {
  #  print.p <- bquote(italic(p)-values<.(signif(2.62, digits=sigdig))%.%10^-14)
  #}
  #text(x= x_control+xoffset, y= control.y,labels = print.p, srt=c(90, 270)[plotnumber])
  
  sample.y <- y_strain[which(order==this.id[1])]
  
  if (sample.y == max(y_strain)) {adjustor <- 0
  } else if (sample.y == min(y_strain)) {adjustor <- 1 
  } else {adjustor <- c(0.4,0.6)[plotnumber]}
  print.p <- formatC(signif(max(p.samples[which(p.samples<0.05)], na.rm = TRUE), digits = sigdig))
  if (as.numeric(print.p)>0) {
    print.p <- strsplit(print.p, split = "e")[[1]]
    if (length(print.p)==1) {
      print.p <- bquote(italic(p)-values<=.(print.p))
    } else {
      print.p <- bquote(italic(p)-values<=.(print.p[1])%.%10^.(as.numeric(print.p[2])))
    }
  } else {
    print.p <- bquote(italic(p)-values<.(signif(2.62, digits=sigdig))%.%10^-14)
  }
  text(x= x_control+xoffset+addme, y= sample.y,labels = print.p, srt=c(90, 270)[plotnumber], adj = adjustor)
  
  
  addme <- addme+0.08
  }
  
  
  #plot strain labels in second plot
  if (plotnumber==1) {
    barchart.oli(means =-1, SE=0.01,labels="", ymin=0, ymax=ymax, xlim=c(1,2), space=1, horiz=TRUE, axes=FALSE)
    for (i in 1:length(get(paste0("means_", thisplot)))) {
      y <- y_strain[i]
      x <- 1.5
      text(x=x, y=y+0.15, labels=labels[i], cex = cex)
      text(x=x, y=y-0.3, labels=paste0("(n=", length(get(paste("PCC", thisplot, order[i], sep = "_"))), ")"), cex=cex)
    }
    mtext(text = "Correlation between fluorescence intensities", side=1, line =3)
    mtext(text = "(Pearson Correlation Coefficient)", side=1, line =4.5)
  }
}



dev.off()
