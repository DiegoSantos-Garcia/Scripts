###Saturation analysis for phylogenetics###
#===============================================================================
#   The idea of this simple test was from Katie Everson, PhD. http://www.kmeverson.org/blog/simple-dna-saturation-plots-in-r
# 
#    Author: Diego Santos-Garcia
#
#    Contact: Contact the author at diego.santos@mail.huji.ac.il or diego.santos@uv.es.
#
#    COPYRIGHT: Copyright (C) 2015  Diego Santos-Garcia.
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(ade4))
suppressPackageStartupMessages(require(seqinr))
suppressPackageStartupMessages(require(optparse))

###Define arguments
option_list = list(
  make_option(c("-i", "--infolder"), action="store", default=".", type='character',
              help="Folder where alingments are stored."),
  make_option(c("-f", "--format"), action="store", default="*fasta.aln$", type='character',
              help='Format pattern present in alingments files. E.g. "*.aln"'),
  make_option(c("-o", "--outfile"), action="store", default="saturation_matrix.tab",
              help="File name for regression coefficient values matrix and saturation plots."),
  make_option(c("-r", "--regression"), action="store", default=0.7,
              help="Select the regression coefficient (0-1) between RAW and corrected distances to filter undesired alingments. Default 0.7.")  
)
opt = parse_args(OptionParser(option_list=option_list))
infolder = opt$i
files_format = opt$f
outfile = opt$o
coeff = opt$r

#Read files and prepare data frame
files = list.files(path=infolder,pattern=files_format)
saturation_matrix = data.frame(matrix(NA,ncol=0,nrow=1))
row.names(saturation_matrix) = c("R2")
pdf(file = paste(outfile,"saturation_plots.pdf",sep="_"))

#Iterates over all the alingment files and calculates the 
for(i in files){
  print(paste(infolder,i,sep =""))
  group = gsub(files_format,"",i)
  print(paste("Working on file ",i,sep=""))
  dat <- read.alignment(file = paste(infolder,i,sep =""),format = "fasta")
  ###Convert to genetic distances###
  dat <- as.DNAbin(dat)
  dist <- dist.dna(dat, model="raw")
  dist.corrected <- dist.dna(dat, model="K80")
  ###Make each plot aln linera regression model
  par(mfrow=c(1,1))
  plot(dist~dist.corrected, pch=20, col="red",xlab="Corrected distance (K80)", ylab="RAW distance", main= i)
  abline(0,1, lty=2)
  abline(lm(dist~dist.corrected), lwd=3)
  lm_coef<-coef(lm(dist~dist.corrected))
  saturation_matrix <- cbind(saturation_matrix,lm_coef[2])
  names(saturation_matrix)[ncol(saturation_matrix)] <- group
  mtext(bquote(r == .(lm_coef[2])))
}
x <- as.data.frame(t(saturation_matrix))
boxplot(x,main="Correlation between RAW and corrected genetic distances",ylab="Correlation (R2)",xlab="Orthologous groups")
satdensity <- density(as.numeric(x$R2))
plot(satdensity,col="blue",main="Correlation values distribution")
x2 <- subset(x,x$R2 >= coeff)
boxplot(x2,main="Correlation between RAW and corrected genetic distances",ylab="Correlation (R2)",xlab="Orthologous groups")
satdensity <- density(as.numeric(x2$R2))
plot(satdensity,col="blue",main="Correlation values distribution")
dev.off()
write.table(x,file=(paste(outfile,"saturation_matrix_RAW.tab",sep="_")),sep="\t")
write.table(x2,file=(paste(outfile,"saturation_matrix_selected.tab",sep="_")),sep="\t")













