#!/usr/bin/env Rscript

# 2013-09-16
# works with R 2.15.2 and ggplot 0.9.3.1
# Check ggplot2 help forums or contact sujai.kumar@gmail.com if something doesn't run
# because of updated programs/packages

#Function to ignore low frequency annotations:

clean.blobs<-function(d,threshold,taxlevel) {
    annotated<-d[d[,taxlevel]!="Not annotated",]
    total<-dim(annotated)[1]
    levels(d[,taxlevel])[which(table(d[,taxlevel])<threshold*total)]<-"Not annotated"
    return(d)
}

#########################################################################

#Load data from file and generate plot:

library(ggplot2)
library(reshape2)

subplotwidth=1000;
subplotheight=1000;

args <- commandArgs(trailingOnly = TRUE)
arg_input_file <- args[1]
arg_ignore_below_prop=as.numeric(args[2])
arg_taxlevel=args[3]

orig <- read.delim(arg_input_file,header=TRUE,sep="\t")
orig <- orig[orig$len>=200,]

# if cov_colnames >1, then create a new column cov_total = sum of all cov columns
if (length(grep("^cov_",colnames(orig))) >1) orig$cov_total = rowSums(orig[,grep("^cov_",colnames(orig))])

cov_colnames=colnames(orig)[grep("^cov_",colnames(orig))]
tax_colnames=colnames(orig)[grep("^taxlevel_",colnames(orig))]

numcols=length(cov_colnames)
 
taxlevel=arg_taxlevel;

m<-melt(orig,id.vars=c("seqid","len","gc",taxlevel),measure.vars=cov_colnames, variable.name="read_set", value.name="cov")

# there aren't many colours available, so to restrict the plot to only 13 colours:
# (thanks to https://github.com/hobrien for the fix)
if (length(levels(m[,taxlevel])) > 14) {
  levels(m[,taxlevel])[which(table(m[,taxlevel])<=sort(as.numeric(table(m[,taxlevel])), decreasing=T)[13])]<-"other"
}

mfilt<-clean.blobs(m,arg_ignore_below_prop,taxlevel)
taxnames=names(sort(table(mfilt[,taxlevel]),decreasing=TRUE))
taxnames=c("Not annotated", taxnames[taxnames != "Not annotated"])
mfilt[,taxlevel] <- factor(mfilt[,taxlevel], levels = taxnames)

png(paste(arg_input_file,taxlevel,"png",sep="."), (numcols * subplotwidth), (1 * subplotheight) + 300, units="px",res=100)

theme_set(theme_bw())
# Paul Tol scheme is well documented at http://www.sron.nl/~pault/colourschemes.pdf - Thank you Paul! Added DDDDDD and 777777 to it
paultol=list(c("#DDDDDD") ,c("#DDDDDD","#4477AA"), c("#DDDDDD","#4477AA","#CC6677"), c("#DDDDDD","#4477AA","#DDCC77","#CC6677"), c("#DDDDDD","#4477AA","#117733","#DDCC77","#CC6677"), c("#DDDDDD","#332288","#88CCEE","#117733","#DDCC77","#CC6677"), c("#DDDDDD","#332288","#88CCEE","#117733","#DDCC77","#CC6677","#AA4499"), c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#DDCC77","#CC6677","#AA4499"), c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#AA4499"), c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499"), c("#DDDDDD","#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#882255","#AA4499"), c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#882255","#AA4499"), c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499"), c("#DDDDDD","#332288","#6699CC","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100","#CC6677","#AA4466","#882255","#AA4499","#777777"))

g<-ggplot() + scale_colour_manual(values=paultol[[length(levels(mfilt[,taxlevel]))]], name="Annotation", limits=levels(mfilt[,taxlevel]) )

for (t in levels(mfilt[,taxlevel])) {
  g <- g + geom_point(data=mfilt[mfilt[,taxlevel]==t,],aes_string(x="gc", y="cov", colour=taxlevel), size=2, alpha=I(1/3))
}

g<-g +
  facet_wrap(~read_set, ncol=numcols) + 
  scale_y_log10(limits=c(5, 10000),breaks = c(1,10,100,1000,10000)) + scale_x_continuous(limits=c(0.2, 0.8),breaks = seq(0,1,.2)) +
  labs(x="GC content", y="Contig abundance") + 
  guides(colour = guide_legend(nrow=3, override.aes = list(alpha = 1,size=10))) + 
  theme (
    strip.text.x = element_text(colour = "black", size = 25, vjust = 0.5),
    axis.text.x  = element_text(colour = "black", size = 25, vjust = 1),
    axis.text.y  = element_text(colour = "black", size = 15, vjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 25, vjust = 0),
    axis.title.y = element_text(colour = "black", size = 25, hjust = 0.5, vjust = 0.5, angle=90),
    legend.text  = element_text(colour = "black", size = 18, vjust = 0),
    legend.title = element_text(colour = "black", size = 20, vjust = 0, hjust = 0, lineheight=1),
    legend.justification=c(1,1), legend.position="bottom", legend.direction="horizontal"
  )
g
dev.off()
