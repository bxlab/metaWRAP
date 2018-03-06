#!/usr/bin/env Rscript

# 2013-09-16
# works with R 2.15.2 and ggplot 0.9.3.1
# Check ggplot2 help forums or contact sujai.kumar@gmail.com if something doesn't run
# because of updated programs/packages

#Function to ignore low frequency annotations:

clean.blobs<-function(d,threshold,taxlevel) {
    annotated<-d[d[,taxlevel]!="Unbinned",]
    total<-dim(annotated)[1]
    levels(d[,taxlevel])[which(table(d[,taxlevel])<threshold*total)]<-"Unbinned"
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
if (length(levels(m[,taxlevel])) > 201) {
  levels(m[,taxlevel])[which(table(m[,taxlevel])<=sort(as.numeric(table(m[,taxlevel])), decreasing=T)[200])]<-"other"
}

mfilt<-clean.blobs(m,arg_ignore_below_prop,taxlevel)
taxnames=names(sort(table(mfilt[,taxlevel]),decreasing=TRUE))
taxnames=c("Unbinned", taxnames[taxnames != "Unbinned"])
mfilt[,taxlevel] <- factor(mfilt[,taxlevel], levels = taxnames)

png(paste(arg_input_file,taxlevel,"png",sep="."), (numcols * subplotwidth), (1 * subplotheight) + 300, units="px",res=100)

theme_set(theme_bw())
random_colors=c("#DDDDDD", "#FA5F70", "#BB24F7", "#EA9143", "#E83EEE", "#E1EE3C", "#70945B", "#CC6B80", "#6B98F3", "#4EB175", "#3EED21", "#916B1A", "#377C4A", "#41BE37", "#F1A9D1", "#71128D", "#01545D", "#35ADC2", "#668A25", "#8BBBF0", "#07E3AE", "#424CDD", "#A7E0A5", "#0A3DC8", "#135826", "#D7304B", "#BE272F", "#ACB103", "#8A1379", "#FB3514", "#999C2A", "#865E66", "#51264C", "#CD05D9", "#57C2F6", "#92B719", "#A98421", "#B40324", "#79B7A0", "#0B3AB7", "#48C323", "#8830D5", "#1E87AA", "#E08E36", "#87AF8F", "#326AFB", "#4C2070", "#ADD98F", "#61FB2E", "#561763", "#DD1E39", "#38392D", "#1BFF29", "#27AC47", "#227A8D", "#AD18C8", "#D19B92", "#DC0577", "#920030", "#ECC6E9", "#BB113A", "#FD76CC", "#7FB5BF", "#EED2B9", "#D6C4DF", "#CB038A", "#3A47AB", "#403F32", "#B75F64", "#7449BA", "#A8BCAF", "#A4DB10", "#7389CE", "#64AA81", "#DE44BB", "#E03223", "#6D234A", "#0124DA", "#FB9E40", "#36221C", "#B42BEF", "#09BCFB", "#65002E", "#5E0AEA", "#5F1274", "#2B7A8C", "#28813E", "#1368B0", "#8DD870", "#F997F9", "#C6C160", "#5FE146", "#BE4647", "#6E195D", "#6C170C", "#93847A", "#4A4E74", "#FD1D18", "#17DFE1", "#4AA396", "#FC3808", "#1F6927", "#AE30C5", "#3F3448", "#78033B", "#98F193", "#38B5EC", "#C2A8BF", "#247C44", "#F4D97E", "#38EF04", "#720E56", "#BC7641", "#DC0064", "#9B5D15", "#FC1E5E", "#DD3717", "#D56288", "#21D33F", "#6DFA34", "#725B55", "#9BAB64", "#DF87C5", "#4B362F", "#262D38", "#B90334", "#3AFF6F", "#0CD849", "#7B781B", "#AB846F", "#B2A7CF", "#070416", "#7D9327", "#ACE44C", "#D2997C", "#9F37EB", "#F34C0F", "#4BB029", "#A1A644", "#BFAB1C", "#5E6607", "#E1E978", "#08E708", "#32FBE1", "#54CF23", "#72750C", "#527ABC", "#88540E", "#0DB5D7", "#4D4278", "#D3DC09", "#E60907", "#F2170D", "#C56456", "#AD71B2", "#B2AE04", "#C3F4ED", "#8185AA", "#576248", "#4C39AA", "#8EA094", "#D22FC5", "#65337A", "#C1421F", "#64C3AC", "#BEEA78", "#B55186", "#203930", "#901115", "#234E7B", "#E58256", "#C5007B", "#BB6654", "#E77D9D", "#01EBC3", "#0D2C2F", "#1E759D", "#3C81D8", "#99283D", "#2BEEC2", "#EDC71A", "#82BA30", "#649623", "#302FD9", "#F1CCBD", "#70BC8A", "#610062", "#5F6845", "#426926", "#34D602", "#A44312", "#5F6895", "#1F7780", "#3B43EF", "#16A94E", "#FE12A7", "#D03C65", "#CA1F26", "#ACD36F", "#7E2048", "#352A79", "#777777")

n <-length(levels(mfilt[,taxlevel]))
colors_new <- random_colors[c(1:n)]

g<-ggplot() + scale_colour_manual(values=colors_new, name="Bin annotation", limits=levels(mfilt[,taxlevel]) )

for (t in levels(mfilt[,taxlevel])) {
  g <- g + geom_point(data=mfilt[mfilt[,taxlevel]==t,],aes_string(x="gc", y="cov", colour=taxlevel), size=2, alpha=I(1/3))
}
#y_axis_breaks = c(1,2,5,10,20,50,100,200,500,1000);

#g<-g+expand_limits(y=c(0,1200))
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
