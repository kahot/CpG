library(colorspace)
library(ggplot2)
library(ggthemes)
library(reshape2)
colors<-qualitative_hcl(6, palette="Dark3")

cpgMF<-read.csv("Data/alldatapoints.csv", stringsAsFactors = F, row.names = 1)
dfA<-cpgMF[,c("Virus","AsynNC_C")]
dfT<-cpgMF[,c("Virus","TsynNC_C")]
colnames(dfA)[2]<-"CpG.Cost.A"
colnames(dfT)[2]<-"CpG.Cost.T"

df<-merge(dfA, dfT, by="Virus")

cpg<-read.csv("Output/Summary/CpG_summary.csv", stringsAsFactors = F, row.names = 1)
cpg<-cpg[,c("Virus", "cg.rho")]

df2<-merge(df, cpg, by="Virus")


plot(log(df2$CpG.Cost.T)~df2$cg.rho)

cor.test(log(df2$CpG.Cost.T), df2$cg.rho,method = "spearman")
#S = 12106, p-value = 0.000503
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#-0.5580438 

cor.test(log(df2$CpG.Cost.A), df2$cg.rho,method = "spearman")
#S = 11238, p-value = 0.006826
#rho = -0.446332 


df3<-melt(df2)
rho<-df2[,c("Virus","cg.rho")]
rho<-rbind(rho,rho)
df.plot<-cbind(df3[1:72,], rho[,2])
colnames(df.plot)[4]<-"CG.rho"
df.plot$variable<-factor(df.plot$variable, levels=c("CpG.Cost.A","CpG.Cost.T" ))

cor.test(df.plot$value, df.plot$CG.rho,method = "spearman" )
#S = 93217, p-value = 8.232e-06
#-0.4987616 

ggplot(df.plot, aes(x=CG.rho, y=value, color=variable))+
    scale_y_continuous(trans="log10",breaks = c(1,2,5, 10,20, 50, 100, 200),labels=c(1,2,5, 10,20, 50, 100, 200))+
    geom_point(stat="identity")+
    theme(legend.title = element_blank())+
    scale_color_manual(values=colors[c(1,5)], labels=c("A","T"))+
    theme_bw()+
    ylab("Cost of CpG creating mutations")+
    xlab("Dinucleotide over/udner-representation")+
    theme(panel.grid.minor=element_blank(),panel.grid.major.x = element_blank())+
    geom_smooth(method='lm',formula=y~x, se=F, size=.3, linetype=2)+
    theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(linetype = 0, size=2)))+
    annotate(geom="text",x=0.75, y=200, label="A: rho=-0.446, P=0.0068", hjust=0, size=2)+
    annotate(geom="text",x=0.75, y=150, label="T: rho=-0.558, P=0.0005", hjust=0, size=2)
ggsave("Output/CPG_correlation.plot1.pdf", width = 5,height = 4)  
     
ggplot(df.plot, aes(x=CG.rho, y=value, color=variable))+
    scale_y_continuous(trans="log10",breaks = c(1,2,5, 10,20, 50, 100, 200),labels=c(1,2,5, 10,20, 50, 100, 200))+
    geom_point(stat="identity")+
    theme(legend.title = element_blank())+
    scale_color_manual(values=colors[c(1,5)], labels=c("A","T"))+
    theme_bw()+
    ylab("Cost of CpG creating mutations")+
    xlab("Dinucleotide over/udner-representation")+
    theme(panel.grid.minor=element_blank())+
    geom_smooth(method='lm',formula=y~x, se=F, size=.3, linetype=2)+
    theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(linetype = 0, size=2)))+
    annotate(geom="text",x=0.75, y=200, label="A: rho=-0.446, P=0.0068", hjust=0, size=2)+
    annotate(geom="text",x=0.75, y=150, label="T: rho=-0.558, P=0.0005", hjust=0, size=2)
ggsave("Output/CPG_correlation.plot2.pdf", width = 5,height = 4)  

ggplot(df.plot, aes(x=CG.rho, y=value, color=variable))+
    scale_y_continuous(trans="log10",breaks = c(1,2,5, 10,20, 50, 100, 200),labels=c(1,2,5, 10,20, 50, 100, 200))+
    geom_point(stat="identity")+
    theme(legend.title = element_blank())+
    scale_color_manual(values=colors[c(1,5)], labels=c("A","T"))+
    theme_bw()+
    ylab("Cost of CpG creating mutations")+
    xlab("Dinucleotide over/udner-representation")+
    theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank())+
    geom_smooth(method='lm',formula=y~x, se=F, size=.3, linetype=2)+
    theme(legend.title = element_blank())+
    guides(color = guide_legend(override.aes = list(linetype = 0, size=2)))+
    annotate(geom="text",x=0.75, y=200, label="A: rho=-0.446, P=0.0068", hjust=0, size=2)+
    annotate(geom="text",x=0.75, y=150, label="T: rho=-0.558, P=0.0005", hjust=0, size=2)
ggsave("Output/CPG_correlation.plot3.pdf", width = 5,height = 4)  
