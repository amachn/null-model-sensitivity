# Supplementary Figure 3

a2 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/real.summary.full.v3.M1N15.txt",head=T, comment.char="") #r01
a4 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/real.summary.full.v3.M1N15.1e-5.txt",head=T, comment.char="") #r01


ids2 = which(a2[,4]<0.000257 & a2[,3]>1)
ids4 = which(a4[,4]<0.000257 & a4[,3]>1)

### r001
pdf("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/SF3_newthresh.pdf",12,6)
par(mfrow=c(1,2),mar=c(5, 5, 4, 2) + 0.1)
	plot(-log10(a2$Pvalue), -log10(a4$Pvalue), main=paste("Pearson correlation ",round(cor(-log10(a2$Pvalue), -log10(a4$Pvalue)),2), sep=""), xlab=expression(paste("-",log[10],"P, T<",10^-8, sep=" ")), ylab=expression(paste("-",log[10],"P, T<",10^-5, sep=" ")),pch=21,bg="skyblue3" ); abline(h=-log10(0.00025),col="tomato",lwd=2); abline(v=-log10(0.00025),col="tomato",lwd=2); text(2,2,length(which(a2$Pvalue>=0.00025 & a4$Pvalue>=0.00025)),col="white", cex=2); text(12,17,length(which(a2$Pvalue<0.00025 & a4$Pvalue<0.00025)),col=1,cex=2); text(12,2,length(which(a2$Pvalue<0.00025 & a4$Pvalue>=0.00025)),col=1,cex=2); text(2,17,length(which(a2$Pvalue>=0.00025 & a4$Pvalue<0.00025)),col=1, cex=2)
	plot((a2$OR[which(a2$Pvalue<0.00025 & a4$Pvalue<0.00025)]), (a4$OR[which(a2$Pvalue<0.00025 & a4$Pvalue<0.00025)]), main=paste("Pearson correlation ",round(cor((a2$OR[which(a2$Pvalue<0.00025 & a4$Pvalue<0.00025)]), (a4$OR[which(a2$Pvalue<0.00025 & a4$Pvalue<0.00025)])),2), sep=""), xlab=expression(paste("OR, T<",10^-8, sep=" ")), ylab=expression(paste("OR, T<",10^-5, sep=" ")), pch=21,bg="skyblue3");
dev.off()


###. Figure 4
library(RColorBrewer)
library(reshape2)
library(ggplot2)


a10 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GOSHIFTER/output_gwas_v2/real.summary.full.h3k27ac.txt22",head=T)
a20 = read.table("/lustre/scratch114/teams/soranzo/users/vi1/gregor/output_gwas_v2_h3k27ac/real.summary.full.txt2",head=T)
a30 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/FGWAS/output_gwas_v2_h3k27ac_NEW/FGWAS.all.summary.h3k27ac.new.txt.2",head=T)
a30$p=1-pchisq(2*(a30[,3]-a30[,7]),1)
a40 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/real.summary.null.annotations.h3k27ac.1e-08.M1N15.txt",head=T, comment.char="")
a50 = read.table("/lustre/scratch114/teams/soranzo/users/vi1/ldsr/output_gwas_h3k27ac/real.summary.full.withbaseline.txt2",head=T, comment.char="")


link = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/histone_modifications/newdata/H3K27ac_garfield_narrowPeak/link_file.txt",head=T)
a20$name = as.character(a20$Bed_File)

fn <- function(x){x[length(x)]}
a10$name =  sapply(strsplit(gsub(".gz","",as.character(a10$names)),split="/",fixed=T),fn)
a10 = a10[grep("narrowPeak",a10$name),]
a30$Annotation = link$Annotation[match(as.numeric(as.character(a30$annotation)) -2, link$Index)] ## checked/ changed index to -2!
a30$Annotation = sapply(strsplit(gsub(".gz","",as.character(a30$Annotation)),split="/",fixed=T),fn)
a40$Annotation = sapply(strsplit(gsub(".gz","",as.character(a40$Annotation)),split="/",fixed=T),fn)
a50$name =  sapply(strsplit(gsub(".gz","",as.character(a50$name)),split="/",fixed=T),fn)
a20$EN=a20[,2]/a20[,3]

dat = data.frame(Enrichment=NULL, Annotation=NULL, P = NULL)
dat = rbind(dat, cbind(unname(a10[,c(5,6,2,1)]),Method="GOSHIFTER"))
dat = rbind(dat, cbind(unname(a20[,c(7,6,4,5)]),Method="GREGOR"))
dat = rbind(dat, cbind(unname(a30[,c(5,9,8,1)]),Method="FGWAS"))
dat = rbind(dat, cbind(unname(a40[,c(3,14,4,19)]),Method="GARFIELD"))
tmp=cbind(unname(a50[,c(10,13,10,11)]),Method="LDSR")
tmp[,3] = 2*(1-pnorm(abs(tmp[,3])))
dat = rbind(dat, tmp)

names(dat) = c("Enrichment","Annotation", "P", "Trait", "Method")
dat$P[which(dat$P==0 & dat$Method=="GOSHIFTER")]=1e-4
dat$neglog10P = -log10(dat$P)

dat2 = dat
dat2$Significance = NA
dat2$Significance[which(dat2$P<0.05)]="0.0005<=P<0.05"
dat2$Significance[which(dat2$P>=0.05)]="P>=0.05"
dat2$Significance[which(dat2$P<=0.00047)]="P<0.0005"
dat2$Significance = factor(dat2$Significance, levels=c("P>=0.05","0.0005<=P<0.05", "P<0.0005"))
dat2$Direction="+"
dat2$Direction[which(dat2$Enrichment<1 & (dat2$Method %in%c("GARFIELD","FGWAS","LDSR")) )]="-"
dat2$Direction = factor(dat2$Direction,levels=c("+","-"))
dat2$neglog10P[which(dat2$Enrichment<1 & dat2$Method!="GOSHIFTER" & dat2$Method!="GREGOR")] = 0
dat2$Index = paste(as.character(dat2$Annotation),as.character(dat2$Trait),sep="-")
dat2$Index2 = paste(dat2$Index,":",dat2$Method,sep="")

dat3 = dat2[which(dat2$Method %in% c("GARFIELD","FGWAS","GOSHIFTER","GREGOR","LDSR") & dat2$Trait %in% levels(dat2$Trait)[-c(1,4,21:24)]),]
dat3 = dat3[which(dat3$Significance=="P<0.0005"),]

agr0 = aggregate(dat3$Trait, list(dat3$Trait,dat3$Annotation),length)
dat3$Numbers2 = agr0$x[match(dat3$Index,paste(as.character(agr0[,2]),as.character(agr0[,1]),sep="-"))]
info = table(dat3$Index)
dat3$Numbers = as.numeric(info[match(dat3$Index,names(info))])
agr1 = aggregate(dat3$Numbers, list(dat3$Method, dat3$Numbers), length)
agr2 = aggregate(dat3$Numbers, list(dat3$Numbers), length)
names(agr1)= c("Method","Numbers","All")
names(agr2)= c("Numbers","All")

agr1$Norm = agr2$All[match(agr1$Numbers, agr2$Numbers)]
agr1$Prop = agr1$All/agr1$Norm
agr1$Method = factor(agr1$Method, levels=c("GREGOR","GARFIELD", "LDSR","FGWAS","GOSHIFTER"))

tis = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/plots/table_mapping_tissues_roadmap.txt",head=T,sep="\t")
tis$Ann = paste(as.character(tis[,1]),"-H3K27ac.narrowPeak",sep="")
dat2$Tissue = tis$Class3[match(dat2$Annotation,tis$Ann)]
range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

tbl0 = NULL
for (trait in levels(dat2$Trait)){
	for (method in levels(dat2$Method)){
		# scale Enrichment values
		dat0 = dat2[which(dat2$Trait==trait & dat2$Method==method), ]
		dat0$Enrichment2 = range01((dat0$Enrichment))

		for (tissue in levels(dat0$Tissue)){
			dat00 = dat0[which(dat0$Trait==trait & dat0$Tissue==tissue & dat0$Method==method), ]
 			tbl0 = rbind(tbl0,dat00[which.min(dat00$P),])
		}	
		print(c(trait,method))
	}
}

levels(tbl0$Trait) = c("T2D","CD","IBD","UC","BMI","HGT","WHR", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","2hrG","FI","HOMA-B","HOMA-IR")

tbl0$Method2 = tbl0$Method
tbl0$Method2[which(tbl0$Significance!="P<0.0005")] = NA
tbl0 = tbl0[which(tbl0$Trait %in% c("BMI","HGT", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC")),]
tbl0$Trait = factor(tbl0$Trait,levels=c("HGT","BMI", "HDL","TG","LDL","TC","MCH","MCV","MCHC","RBC","HGB","PCV","PLT","MPV","CD","IBD","UC","FPI","FG","HbA1C","SCZ"))
tbl0$Method2 = factor(tbl0$Method2,levels=c("GARFIELD","GREGOR","LDSR","FGWAS","GOSHIFTER"))
tbl0$Method = factor(tbl0$Method,levels=c("GARFIELD","GREGOR","LDSR","FGWAS","GOSHIFTER"))
tbl0$Tissue = gsub("_"," ",tbl0$Tissue)


pdf("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4c_newthresh2.pdf",20,6)
p<-ggplot(tbl0) +geom_tile(aes(Method,Tissue),fill="grey90")+geom_point(aes(Method,Tissue,colour=Method2, size=Enrichment2), shape=15)+ theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks = element_blank(),axis.text.y=element_text(size=13), axis.title=element_text(size=18),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=11),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 6,colour="white"))+ylab("Tissue")+scale_colour_manual(values=c("red","orange","lightgoldenrod2","seagreen3","steelblue3"),na.value="grey70")+guides(fill=FALSE) +facet_grid(Tissue~Trait,scales="free_y")+scale_size_area(max_size=4)
print(p)
dev.off()

system("evince /lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4c.pdf &")


#### extract supplementary table of results

tbl0f = dat2[,c(2,11,4,6,5)]
levels(tbl0f$Trait) = c("T2D","CD","IBD","UC","BMI","HGT","WHR", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","2hrG","FI","HOMA-B","HOMA-IR")
tbl0f = tbl0f[which(tbl0f$Trait %in% c("BMI","HGT", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC")),]
tbl0f$Tissue = gsub("_"," ",tbl0f$Tissue)

library(reshape2)
data_wide <- dcast(tbl0f, Annotation+Tissue+Trait ~ Method, value.var="neglog10P")

write.table(data_wide,file="/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/ST5b.txt", sep="\t",row.names=FALSE,col.names=TRUE, quote=FALSE, append=FALSE)

####



a10 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GOSHIFTER/output_gwas_v2/real.summary.full.h3k4me3.txt22",head=T)
a20 = read.table("/lustre/scratch114/teams/soranzo/users/vi1/gregor/output_gwas_v2_h3k4me3/real.summary.full.txt2",head=T)
a30 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/FGWAS/output_gwas_v2_h3k4me3_NEW/FGWAS.all.summary.h3k4me3.new.txt.2",head=T)
a40 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/real.summary.null.annotations.h3k4me3.1e-08.M1N15.txt",head=T, comment.char="") 
a50 = read.table("/lustre/scratch114/teams/soranzo/users/vi1/ldsr/output_gwas_h3k4me3/real.summary.full.withbaseline.txt2",head=T, comment.char="")
a30$p=1-pchisq(2*(a30[,3]-a30[,7]),1)

## unify annotation names
link = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/histone_modifications/newdata/H3K4me3_garfield_narrowPeak/link_file.txt",head=T)
a20$name = as.character(a20$Bed_File)

fn <- function(x){x[length(x)]}
a10$name =  sapply(strsplit(gsub(".gz","",as.character(a10$names)),split="/",fixed=T),fn)
a10 = a10[grep("narrowPeak",a10$name),]
a30$Annotation = link$Annotation[match(as.numeric(as.character(a30$annotation)) -2, link$Index)] ## checked/ changed index to -2!
a30$Annotation = sapply(strsplit(gsub(".gz","",as.character(a30$Annotation)),split="/",fixed=T),fn)
a40$Annotation = sapply(strsplit(gsub(".gz","",as.character(a40$Annotation)),split="/",fixed=T),fn)
a50$name =  sapply(strsplit(gsub(".gz","",as.character(a50$name)),split="/",fixed=T),fn)
a20$EN=a20[,2]/a20[,3]

dat = data.frame(Enrichment=NULL, Annotation=NULL, P = NULL)
dat = rbind(dat, cbind(unname(a10[,c(5,6,2,1)]),Method="GOSHIFTER"))
dat = rbind(dat, cbind(unname(a20[,c(7,6,4,5)]),Method="GREGOR"))
dat = rbind(dat, cbind(unname(a30[,c(5,9,8,1)]),Method="FGWAS"))
dat = rbind(dat, cbind(unname(a40[,c(3,14,4,19)]),Method="GARFIELD"))

tmp=cbind(unname(a50[,c(10,13,10,11)]),Method="LDSR")
tmp[,3] = 2*(1-pnorm(abs(tmp[,3])))
dat = rbind(dat, tmp)


names(dat) = c("Enrichment","Annotation", "P", "Trait", "Method")
dat$P[which(dat$P==0 & dat$Method=="GOSHIFTER")]=1e-4
dat$neglog10P = -log10(dat$P)

dat2 = dat

dat2$Significance = NA
dat2$Significance[which(dat2$P<0.05)]="0.0005<=P<0.05"
dat2$Significance[which(dat2$P>=0.05)]="P>=0.05"
dat2$Significance[which(dat2$P<=0.00047)]="P<0.0005"
dat2$Significance = factor(dat2$Significance, levels=c("P>=0.05","0.0005<=P<0.05", "P<0.0005"))
dat2$Direction="+"
dat2$Direction[which(dat2$Enrichment<1 & (dat2$Method %in%c("GARFIELD","FGWAS","LDSR")) )]="-"
dat2$Direction = factor(dat2$Direction,levels=c("+","-"))

dat2$neglog10P[which(dat2$Enrichment<1 & dat2$Method!="GOSHIFTER" & dat2$Method!="GREGOR")] = 0 #### Maybe comes from here????

library(RColorBrewer)
library(reshape2)
library(ggplot2)

dat2$Index = paste(as.character(dat2$Annotation),as.character(dat2$Trait),sep="-")
dat2$Index2 = paste(dat2$Index,":",dat2$Method,sep="")

dat3 = dat2[which(dat2$Method %in% c("GARFIELD","FGWAS","GOSHIFTER","GREGOR","LDSR") & dat2$Trait %in% levels(dat2$Trait)[-c(1,4,21:24)]),]
dat3 = dat3[which(dat3$Significance=="P<0.0005"),]

agr0 = aggregate(dat3$Trait, list(dat3$Trait,dat3$Annotation),length)
dat3$Numbers2 = agr0$x[match(dat3$Index,paste(as.character(agr0[,2]),as.character(agr0[,1]),sep="-"))]
info = table(dat3$Index)
dat3$Numbers = as.numeric(info[match(dat3$Index,names(info))])
agr1 = aggregate(dat3$Numbers, list(dat3$Method, dat3$Numbers), length)
agr2 = aggregate(dat3$Numbers, list(dat3$Numbers), length)
names(agr1)= c("Method","Numbers","All")
names(agr2)= c("Numbers","All")
agr1$Norm = agr2$All[match(agr1$Numbers, agr2$Numbers)]
agr1$Prop = agr1$All/agr1$Norm
agr1$Method = factor(agr1$Method, levels=c("GREGOR","GARFIELD", "LDSR","FGWAS","GOSHIFTER"))

tis = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/plots/table_mapping_tissues_roadmap.txt",head=T,sep="\t")
tis$Ann = paste(as.character(tis[,1]),"-H3K4me3.narrowPeak",sep="")
dat2$Tissue = tis$Class3[match(dat2$Annotation,tis$Ann)]
range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

tbl0 = NULL
for (trait in levels(dat2$Trait)){
	for (method in levels(dat2$Method)){
		# scale Enrichment values
		dat0 = dat2[which(dat2$Trait==trait & dat2$Method==method), ]
		dat0$Enrichment2 = range01((dat0$Enrichment))

		for (tissue in levels(dat0$Tissue)){
			dat00 = dat0[which(dat0$Trait==trait & dat0$Tissue==tissue & dat0$Method==method), ]
			tbl0 = rbind(tbl0,dat00[which.min(dat00$P),])
		}
		print(c(trait,method))
	}
}
levels(tbl0$Trait) = c("CDold","T2D","CD","IBD","UC","BMI","HGT","WHR", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","UCold","2hrG","FI","HOMA-B","HOMA-IR")

tbl0$Method2 = tbl0$Method
tbl0$Method2[which(tbl0$Significance!="P<0.0005")] = NA
tbl0 = tbl0[which(tbl0$Trait %in% c("BMI","HGT", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC")),]
tbl0$Trait = factor(tbl0$Trait,levels=c("HGT","BMI", "HDL","TG","LDL","TC","MCH","MCV","MCHC","RBC","HGB","PCV","PLT","MPV","CD","IBD","UC","FPI","FG","HbA1C","SCZ"))
tbl0$Method2 = factor(tbl0$Method2,levels=c("GARFIELD","GREGOR","LDSR","FGWAS","GOSHIFTER"))
tbl0$Method = factor(tbl0$Method,levels=c("GARFIELD","GREGOR","LDSR","FGWAS","GOSHIFTER"))
tbl0$Tissue = gsub("_"," ",tbl0$Tissue)

pdf("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4d_newthresh2.pdf",20,6)
p<-ggplot(tbl0) +geom_tile(aes(Method,Tissue),fill="grey90")+geom_point(aes(Method,Tissue,colour=Method2, size=Enrichment2), shape=15)+ theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks = element_blank(),axis.text.y=element_text(size=13), axis.title=element_text(size=18),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=11),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 6,colour="white"))+ylab("Tissue")+scale_colour_manual(values=c("red","orange","lightgoldenrod2","seagreen3","steelblue3"),na.value="grey70")+guides(fill=FALSE) +facet_grid(Tissue~Trait,scales="free_y")+scale_size_area(max_size=4)
print(p)
dev.off()

system("evince /lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4d.pdf &")

tbl0f = dat2[,c(2,11,4,6,5)]
levels(tbl0f$Trait) = c("CDold","T2D","CD","IBD","UC","BMI","HGT","WHR", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","UCold","2hrG","FI","HOMA-B","HOMA-IR")
tbl0f = tbl0f[which(tbl0f$Trait %in% c("BMI","HGT", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC")),]
tbl0f$Tissue = gsub("_"," ",tbl0f$Tissue)

library(reshape2)
data_wide <- dcast(tbl0f, Annotation+Tissue+Trait ~ Method, value.var="neglog10P")

write.table(data_wide,file="/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/ST5c.txt", sep="\t",row.names=FALSE,col.names=TRUE, quote=FALSE, append=FALSE)


########
# Fig 4 a/b

a10 = rbind(read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GOSHIFTER/output_gwas_v2/real.summary.full.txt",head=T),
read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GOSHIFTER/output_gwas_v2/real.summary.full.part2.txt",head=T))
a20 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GREGOR/output_gwas_v2/real.summary.full.txt2",head=T)
a30 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/FGWAS/output_gwas/FGWAS.all.summary.NEW.txt.2",head=T)
a40 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/real.summary.full.v3.M1N15.txt",head=T, comment.char="") #r01
a50 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/LDSR/output_gwas_v2/real.summary.full.withbaseline.v2.txt22",head=T, comment.char="")
a405 = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/real.summary.full.v3.M1N15.1e-5.txt",head=T, comment.char="") #r01
a30$p=1-pchisq(2*(a30[,3]-a30[,7]),1)

## unify annotation names
link = read.table("/lustre/scratch119/humgen/projects/uk10k/users/vi1/garfield-data0/annotation/link_file.txt",head=T)
link = link[which(link$Category=="Hotspots"),]
a20$name = as.character(a20$Bed_File)
fn <- function(x){x[length(x)]}
a10$name =  sapply(strsplit(gsub(".gz","",as.character(a10$names)),split="/",fixed=T),fn)
a30$Annotation = link$Annotation[match(as.numeric(as.character(a30$annotation)) -1, link$Index)] ## checked/ changed index to -2!
a40$Annotation = as.character(a40$Annotation)
a405$Annotation = as.character(a405$Annotation)
a50$name =  sapply(strsplit(gsub(".gz","",as.character(a50$name)),split="/",fixed=T),fn)
a20$EN=a20[,2]/a20[,3]

a30$Annotation = paste(as.character(a30$Annotation),".0",sep="")
a40$Annotation = paste(as.character(a40$Annotation),".0",sep="")
a405$Annotation = paste(as.character(a405$Annotation),".0",sep="")

dat = data.frame(Enrichment=NULL, Annotation=NULL, P = NULL)
dat = rbind(dat, cbind(unname(a10[,c(5,6,2,1)]),Method="GOSHIFTER"))
dat = rbind(dat, cbind(unname(a20[,c(7,6,4,5)]),Method="GREGOR"))
dat = rbind(dat, cbind(unname(a30[,c(5,9,8,1)]),Method="FGWAS"))
dat = rbind(dat, cbind(unname(a40[,c(3,14,4,19)]),Method="GARFIELD"))
tmp=cbind(unname(a50[,c(10,13,10,11)]),Method="LDSR")
tmp[,3] = 2*(1-pnorm(abs(tmp[,3])))
dat = rbind(dat, tmp)

names(dat) = c("Enrichment","Annotation", "P", "Trait", "Method")
dat$P[which(dat$P==0 & dat$Method=="GOSHIFTER")]=1e-4
dat$neglog10P = -log10(dat$P)

dat2 = dat
dat2$Significance = NA
dat2$Significance[which(dat2$P<0.05)]="0.00025<=P<0.05"
dat2$Significance[which(dat2$P>=0.05)]="P>=0.05"
dat2$Significance[which(dat2$P<=0.000257)]="P<0.00025"
dat2$Significance = factor(dat2$Significance, levels=c("P>=0.05","0.00025<=P<0.05", "P<0.00025"))
dat2$Direction="+"
dat2$Direction[which(dat2$Enrichment<1 & (dat2$Method %in%c("GARFIELD","FGWAS","LDSR")) )]="-"
dat2$Direction = factor(dat2$Direction,levels=c("+","-"))
dat2$neglog10P[which(dat2$Enrichment<1 & dat2$Method!="GOSHIFTER" & dat2$Method!="GREGOR")] = 0

library(RColorBrewer)
library(reshape2)
library(ggplot2)

dat2$Index = paste(as.character(dat2$Annotation),as.character(dat2$Trait),sep="-")
dat2$Index2 = paste(dat2$Index,":",dat2$Method,sep="")

dat3 = dat2[which(dat2$Method %in% c("GARFIELD","FGWAS","GOSHIFTER","GREGOR","LDSR") & dat2$Trait %in% levels(dat2$Trait)[-c(1,4,24:29)]),]
dat3 = dat3[which(dat3$Significance=="P<0.00025"),]

agr0 = aggregate(dat3$Trait, list(dat3$Trait,dat3$Annotation),length)
dat3$Numbers2 = agr0$x[match(dat3$Index,paste(as.character(agr0[,2]),as.character(agr0[,1]),sep="-"))]
info = table(dat3$Index)
dat3$Numbers = as.numeric(info[match(dat3$Index,names(info))])
agr1 = aggregate(dat3$Numbers, list(dat3$Method, dat3$Numbers), length)
agr2 = aggregate(dat3$Numbers, list(dat3$Numbers), length)
names(agr1)= c("Method","Numbers","All")
names(agr2)= c("Numbers","All")
agr1$Norm = agr2$All[match(agr1$Numbers, agr2$Numbers)]
agr1$Prop = agr1$All/agr1$Norm
agr1$Method = factor(agr1$Method, levels=c("GREGOR","GARFIELD", "LDSR","FGWAS","GOSHIFTER"))

pdf("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4a_newthresh.pdf", 12,5)
ggplot(agr1,aes(Numbers,Prop*Numbers, fill=Method))+geom_bar(stat="identity",position="dodge",colour="black")+ theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=13), axis.title=element_text(size=14,face="bold"))+ylab("Proportion of enrichments attributable to method")+xlab("Number of methods declaring enichment")+scale_fill_manual(values=c("orange","red","lightgoldenrod2","lightgreen","steelblue3"))+ theme(legend.justification=c(0,1), legend.position=c(0.0,0.9))+guides(colour=FALSE)
dev.off()


####
dat2$Tissue = link[match(dat2$Annotation,paste(link$Annotation,".0",sep="")),4]
range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

tbl0 = NULL
for (trait in levels(dat2$Trait)){
	for (method in levels(dat2$Method)){
		# scale Enrichment values
		dat0 = dat2[which(dat2$Trait==trait & dat2$Method==method), ]
		dat0$Enrichment2 = range01((dat0$Enrichment))
		
		for (tissue in levels(dat0$Tissue)){
			dat00 = dat0[which(dat0$Trait==trait & dat0$Tissue==tissue & dat0$Method==method), ]
			tbl0 = rbind(tbl0,dat00[which.min(dat00$P),])
		}
		print(c(trait,method))
	}
}
levels(tbl0$Trait) = c("T2D","BMI","HGT","WHR", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC","2hrG","FI","HOMA-B","HOMA-IR","DBP","SBP")

tbl0$Method2 = tbl0$Method
tbl0$Method2[which(tbl0$Significance!="P<0.00025")] = NA
tbl0 = tbl0[which(tbl0$Trait %in% c("BMI","HGT", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC")),]
tbl0$Trait = factor(tbl0$Trait,levels=c("HGT","BMI", "HDL","TG","LDL","TC","MCH","MCV","MCHC","RBC","HGB","PCV","PLT","MPV","CD","IBD","UC","FPI","FG","HbA1C","SCZ"))
tbl0$Method2 = factor(tbl0$Method2,levels=c("GARFIELD","GREGOR","LDSR","FGWAS","GOSHIFTER"))
tbl0$Method = factor(tbl0$Method,levels=c("GARFIELD","GREGOR","LDSR","FGWAS","GOSHIFTER"))
tbl0$Tissue = gsub("_"," ",tbl0$Tissue)

pdf("/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4b_newthresh.pdf",20,12)
ggplot(tbl0) +geom_tile(aes(Method,Tissue),fill="grey90")+geom_point(aes(Method,Tissue,colour=Method2, size=Enrichment2), shape=15)+ theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks = element_blank(),axis.text.y=element_text(size=13), axis.title=element_text(size=18),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=11),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 6,colour="white"))+ylab("Tissue")+scale_colour_manual(values=c("red","orange","lightgoldenrod2","seagreen3","steelblue3"),na.value="grey70")+guides(fill=FALSE) +facet_grid(Tissue~Trait,scales="free")+scale_size_area(max_size=4)
dev.off()
system("evince /lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/F4b.pdf &")


#### extract supplementary table of results

tbl0f = dat2[,c(2,11,4,6,5)]
levels(tbl0f$Trait) = c("T2D","BMI","HGT","WHR", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC","2hrG","FI","HOMA-B","HOMA-IR","DBP","SBP")
tbl0f = tbl0f[which(tbl0f$Trait %in% c("BMI","HGT", "HDL","HGB","LDL","FG","HbA1C","FPI","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","CD","IBD","UC")),]
tbl0f$Tissue = gsub("_"," ",tbl0f$Tissue)

library(reshape2)
data_wide <- dcast(tbl0f, Annotation+Tissue+Trait ~ Method, value.var="neglog10P")

write.table(data_wide,file="/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/ST5a_newthresh.txt", sep="\t",row.names=FALSE,col.names=TRUE, quote=FALSE, append=FALSE)

aa = aggregate(tbl0f$neglog10P>(-log10(0.000257)), list(tbl0f$Method, tbl0f$Trait), sum)

 aggregate(aa$x,list(aa$Group.1), summary)
    Group.1   x.Min. x.1st Qu. x.Median   x.Mean x.3rd Qu.   x.Max.
1 GOSHIFTER   0.0000    0.0000   0.0000   0.5238    0.0000   5.0000
2    GREGOR   0.0000    6.0000  24.0000  83.4300  136.0000 398.0000
3     FGWAS   0.0000    0.0000   5.0000  53.8600   51.0000 327.0000
4  GARFIELD   0.0000    0.0000  10.0000  56.2400   88.0000 364.0000
5      LDSR   0.0000    0.0000   5.0000  17.3300   20.0000 144.0000


##############################

a = read.table("../output_v2/real.summary.real.annotations.full.cap.r001.m1,n15,t5.ST4.txt", head=T)

levels(a$trait) = c("CDold","T2D","CD","IBD","UC","BMI","HGT","WHR", "HDL","HGB","DBP","SBP","LDL","2hrG","FG","HbA1C","FI","FPI","HOMA-B","HOMA-IR","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG","UCold")
a = a[which(a$trait %in% c("T2D","CD","IBD","UC","BMI","HGT","WHR", "HDL","HGB","DBP","SBP","LDL","2hrG","FG","HbA1C","FI","FPI","HOMA-B","HOMA-IR","MCH","MCHC","MCV","MPV","PCV","PLT","RBC","SCZ","TC","TG")),]
a$trait = as.factor(as.character(a$trait))

pt=0.000257

tbl2f = a[which(a$Pvalue<pt),c(1,14,15,16,19,2,3,4,5,6,7,8)]
names(tbl2f)[1:5]=c("Index","Annotation", "Celltype","Tissue", "Trait")
write.table(tbl2f,file = "/lustre/scratch119/humgen/projects/uk10k/users/vi1/enrichment_method_comparison/GARFIELD/output_v2/figures/ST4_newthresh.txt", col.names=TRUE,row.names=FALSE, append=FALSE,quote=FALSE, sep="\t")


