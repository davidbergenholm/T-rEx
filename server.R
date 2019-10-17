#
# This is a Shiny web application called T-rex. 
#
# Authors: David Bergenholm, Christoph Börlin, Petter Holland & Jens Nielsen
# Chalmers University of Technology
# Department of Biology and Bioengineering

#Set up the colors, conditions and the TFs used
#For each new TF added a new color must be set and added in rinputs
myColors.TF<-c("#a6cee3","#1f78b4","#ff9c63","#f16913","#b2df8a","#ffe60c","#33a02c","#3ab4f4","#5200aa","#e31a1c","#fdbf6f","#ff7f00","#115e0c","#cab2d6","#7fcdbb","#6a3d9a","#b0b227","#e20002","#dd3497","#b15928")
name.TF<-c("Cat8","Cbf1","Ert1","Gcn4","Gcr1","Gcr2","Hap1","Ino2","Ino4","Leu3","Oaf1","Pip2","Rds2","Rgt1","Rtg1","Rtg3","Sip4","Stb5","Sut1","Tye7")
name.Cond<-c("Glu","Nit","Eth","Ana")
colorset.TF<-as.character(t(myColors.TF))
names(colorset.TF)<-as.character(t(name.TF))
#Make TF list
TF_list<-c()
for(tf in name.TF){
  TF_list[tf]<-list(tf=tf)
}

#Load TSS start position
TSS.start<-read.csv(paste("TF_data_files/Resources/","190704_TSSData.tsv",sep=""), sep="\t", header=TRUE)

#Get the current genelist
geneList<-data.frame(TSS.start$GeneName, TSS.start$Gene)
colnames(geneList)<-c("Gene","GeneC")
gene.Systematic<-data.frame(geneList[,1])
colnames(gene.Systematic)<-c("GeneS")
gene.Common<-data.frame(geneList[,2])
colnames(gene.Common)<-c("GeneC")

#Load metabolic genes from Yeast8
metabolicgenes <-tryCatch({read.csv(paste("TF_data_files/Resources/","Yeast8_genes.csv",sep=""),sep=";", header=TRUE)},
                    error=function(file){
                      dummyvec<-data.frame(t(rep(NA, 1)))
                      colnames(dummyvec)<-c("Gene")
                      return(dummyvec) })

#Load TATA position, realigned from S288C to Cen.PK
TATA.file <-tryCatch({read.cssv(paste("TF_data_files/Resources/","TATA_pos.csv",sep=""),sep=";", header=TRUE)},
                    error=function(file){
                      dummyvec<-data.frame(t(rep(NA, 3)))
                      colnames(dummyvec)<-c("gene_id", "gene_common","TATA_cen_pos")
                      return(dummyvec) })

#Load Sequence for each Gene, 2000 bp centered around TSS
Geneseq <-tryCatch({read.csv(paste("TF_data_files/Resources/","TSS_anotation_seq.bed",sep=""),sep="\t",header = FALSE )},
                    error=function(file){
                      dummyvec<-data.frame(t(rep(NA, 3)))
                      colnames(dummyvec)<-c("Gene","Strand","Seq")
                      return(dummyvec) })

#Load TPM for genes in each condition
TPM.data <-tryCatch({read.csv(paste("TF_data_files/Resources/","TPM.csv",sep=""),sep=";")},
                    error=function(file){
                      NRcond<-t(rep(2, (length(name.Cond))))
                      dummyvec<-data.frame(geneList, NRcond)
                      colnames(dummyvec)<-c("Gene","CommonName",name.Cond)
                      return(dummyvec) })

#create dataframes for TFs
x  <- matrix(1:2000)
x.three<-seq(1,2000,3)

####Load GO-terms located in the gogenes file, if not then load from Ensembl
mtry <- try(read.csv(paste("TF_data_files/Resources/","gogenes.csv",sep=""),sep=";"), 
            silent = TRUE)
if (class(mtry) != "try-error") {
  go_bio<-read.csv(paste("TF_data_files/Resources/","gogenes.csv",sep=""),sep=";")
} else {
  message("File doesn't exist, loading from Ensembl")
  ensembl=useMart("ENSEMBL_MART_ENSEMBL", "scerevisiae_gene_ensembl", host="www.ensembl.org")
  gogenes<-getBM(attributes=c('ensembl_gene_id','name_1006','namespace_1003'),mart=ensembl)
  
  go_bio<-gogenes[grep("biological_process",gogenes$namespace_1003),] # Select only Biological Process
  go_bio<-suppressMessages(remove.vars(go_bio,c("namespace_1003")))
  
}
go_bio_unique<-data.frame(unique(go_bio$name_1006))

#Load GEM peaks 
for (i in name.TF){
  if(i==name.TF[1]){
    file.name<-list.files(path="TF_data_files/Data/",pattern=paste(i,"(.*)_GEManalysis_(.*).csv$",sep=""))
    datafile <-tryCatch({read.csv(paste0("TF_data_files/Data/", file.name),sep=",")},
                           error=function(file){
                             dummyvec<-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
                             colnames(dummyvec)<-c("X","Gene","Chr","Pos","Strand","DistanceTSS","Strength","Condition","TF")
                             return(dummyvec) })
    datafile$TF<-i
  }else{
    file.name<-list.files(path="TF_data_files/Data/",pattern=paste(i,"(.*)_GEManalysis_(.*).csv$",sep=""))
    datafile2 <-tryCatch({read.csv(paste0("TF_data_files/Data/", file.name),sep=",")},
                        error=function(file){
                          dummyvec<-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
                          colnames(dummyvec)<-c("X","Gene","Chr","Pos","Strand","DistanceTSS","Strength","Condition","TF")
                          return(dummyvec) })
    datafile2$TF<-i
    datafile<-rbind(datafile,datafile2)}
    datafile_Peaks<-datafile
  
}

#Generates the Transcription Factor Binding Data
Reads_func<- function(tf,cond,data.New.TF) {
  if(tf=="New_TF"){
    p=data.New.TF[[paste("WigLike",cond,sep="")]]
    return(p)
  }else{
    file.name<-list.files(path="TF_data_files/Data/",pattern=paste(tf,"_",cond,"_ol_combRep_geneAssigned(.*).wigLike$",sep=""))
    p<-tryCatch({read.csv(paste0("TF_data_files/Data/", file.name),sep="\t",header=FALSE)},
                         error=function(file){
                           dummyvec<-data.frame(NA,NA,NA)
                           return(dummyvec) })
  }}

#Generate the motif, sequence map and the peak distribution images
motif_func<-function(input.TF, condition){
  out<-list(
    src = paste("TF_data_files/Data/",input.TF,"_",condition,"_logo1.png",sep=""),
    contentType = "image/png",
    alt = "motif",
    height    = 100,
    units     = "in"
  )
  return(out)
}
seq_map_func<-function(input.TF, condition){
  out<-list(
    src = paste("TF_data_files/Data/",input.TF,"_",condition,"_SeqMap.png",sep=""),
    contentType = "image/png",
    alt = "seq_map",
    height    = 100,
    units     = "in"
  )
  return(out)
}
peak_dist_func<-function(input.TF, condition){
   name.file<-list.files(path="TF_data_Files/Data", pattern=paste(input.TF,"_",condition,"_PeakHistogram(.*).png", sep=""))

  out<-list(
    
    src = paste("TF_data_files/Data/",name.file,sep=""),
    contentType = "image/png",
    alt = "peak_dist",
    height    = 170,
    units     = "in"
  )
  return(out)
}


Stat_plot_func<-function(val1,downloadVal,conditionname,goterm,name.TF,Test,New.TF,nclu){
  GOterm.name<-unlist(strsplit(tolower(goterm), "[+]"))
  GOterm.Data<-goterm_func(GOterm.name)
  totalgenes<-nrow(GOterm.Data)
  ####The data

  Stat.Data<-Stat_data_func(name.TF,conditionname,New.TF$Targets)
  Stat.GOterm.data<-data_func(Stat.Data,conditionname,GOterm.Data,val1)
  Treated.Data<-data_treatment_func(Stat.GOterm.data$data1, conditionname,name.TF)
  txtstr<-paste( "Selected genes", nrow(Treated.Data$x),"of total",totalgenes)
  data<-cbind(Treated.Data$y,Treated.Data$x)
  ###Fisher
  if(Test=="Fisher"){
    p1<-fishers_test_func(Treated.Data$x, txtstr,downloadVal)
                 
  }
  if(Test=="Heatmap"){
    ##### Only data heatmap
    p1<-heatmap_func(data,txtstr,downloadVal)
    
  }
  if(Test=="Network"){
    p1<-net_func(Treated.Data$x, txtstr)
    
  }
  if(Test=="Cluster"){
    p1<-cluster_func(Treated.Data$x,downloadVal, nclu)
    
  }
  if(Test=="Linear Model"){
    # datalin<-cbind(Treated.Data$y,Treated.Data$x)
    #zero interaction
    p1<- model_zero_func(data, txtstr,downloadVal,name.TF)
  }
  
  p1
}
##Find the GO-term and the genes
goterm_func <- function(datain) {
  ####Locates the GO-terms included in the users searchterm
  pathways.name<-datain
  pathways.go <- filter(go_bio, grepl(paste(pathways.name, collapse="|"),name_1006))
  
  pathways.gene<-data.frame(unique(pathways.go$ensembl_gene_id))
  colnames(pathways.gene)<-"Gene"
  return(pathways.gene)
}
#Generate data for statistical analysis 
Stat_data_func<-function(tf,cond,data.New.TF) {
  data_out<-geneList
  for (i in tf){
    if(i=="New_TF"){
      datafile<-data.New.TF
      condTF<-datafile[names(datafile) == cond] 
      datafile.2<-data.frame(datafile$X,condTF)
      colnames(datafile.2)<-c("GeneC",i)
      data_out<-merge(data_out,datafile.2, by="GeneC", all=TRUE)
    }else{
      name<-list.files(path="TF_data_files/Data/",pattern=paste(i,"_geneTargetList_(.*).csv$",sep=""))
      datafile <-tryCatch({read.csv(paste0("TF_data_files/Data/",name),sep=",")},
                         error=function(file){
                           dummyvec<-data.frame(t(rep(NA, (length(name.Cond)+1))))
                           colnames(dummyvec)<-c("X",name.Cond)
                           return(dummyvec) })
      
      condTF<-datafile[names(datafile) == cond] 
      datafile.2<-data.frame(datafile$X,condTF)
      colnames(datafile.2)<-c("GeneC",i)
      data_out<-merge(data_out,datafile.2, by="GeneC", all=TRUE)
    }}
  return(data_out)
  
}
data_func<- function(Stat.Data,datain,pathway.val,val1) {
  ####Generates the data used, if Yeast8 and dubious is or isn't selected
  data1<-Stat.Data
  ####Remove non-metabolic genes
  if (val1==1){
    data1<-merge(data1,metabolicgenes,by="Gene")
  }
  ####Merge with pathway genes
  data1<-merge(data1,pathway.val,by="Gene")
  if(is.null(data1)){data1<-c("Gene"="empty")}
  return(list("data1"=data1))
}
data_treatment_func <- function(data1, condname,name.TF) {

  data1<-merge(data1, TPM.data, by="Gene")
  y1<-data1[condname]
  x1<-data1[,name.TF]
  rownames(x1)<-data1$CommonName
  ####Only include genes which has a TPM above 1
  x1<-x1[y1>1,]
  x1[is.na(x1)]<-0
  y1<-y1[y1>1]
  y1<-log2(y1)
  ####only include genes which has at least one binding 
  y1<-y1[rowSums(x1)>0]
  x1<-x1[rowSums(x1)>0,]
  x1[is.na(x1)]<-0
  return(list("y"=y1,"x"=x1))
}
###Functions for generating Fisherplots, Heatmap, Net plot, Linear Model or Clusters
fishers_test_func <- function(datain, txtstr1,outputval) {
  df1<-as.data.frame(datain)
  df1[df1>0]<-1
  df1[df1<0]<-0
  fishertable1<-data.frame()
  fishertableodds1<-data.frame()
  #### loop over the different TFs
  for (i in 1:length(df1))  {
    for (o in 1:length(df1)) {
     
      first <- ifelse(df1[i]>0,"Bound","NotBound")
      second <- ifelse(df1[o]>0,"Bound","NotBound")
      
      if (any(grepl("^NotBound$", first))==FALSE) {
        rand.not<-cbind("NotBound")
        rownames(rand.not)<-"Rand"
        first<-rbind(first,rand.not)
        second<-rbind(second,rand.not) 
       
      }
      if (any(grepl("^NotBound$", second))==FALSE) {
        rand.not<-cbind("NotBound")
        rownames(rand.not)<-"Rand"
        first<-rbind(first,rand.not)
        second<-rbind(second,rand.not) 
        
      }
      if (any(grepl("^Bound$", first))==TRUE) {
        if (any(grepl("^Bound$", second))==TRUE) {
          fishers<-fisher.test(table(first,second))
          fishertable1[i,o]<-fishers$p.value
          fishertableodds1[i,o]<-fishers$estimate
        } 
        else{
          fishertable1[i,o]<-1
          fishertableodds1[i,o]<-0
        }
      }
      else{
        fishertable1[i,o]<-1
        fishertableodds1[i,o]<-0
      }
    }
  }
  colnames(fishertable1)<-colnames(df1);
  colnames(fishertableodds1)<-colnames(df1);
  rownames(fishertable1)<-colnames(df1);
  rownames(fishertableodds1)<-colnames(df1);
  
  #### The diagonal becomes infinite if the transcription factor correlates with it self, thus we set a dummy OR to make scaling easier
  fishertableodds1[mapply(is.infinite, fishertableodds1)] <- 100
  
  
  #### Create corr matrix and plot corr heatmap
  sign.lev <- 0.001
  OR <- as.matrix(fishertableodds1)
  pval <- as.matrix(fishertable1)
  #### If OR < 1 (co-located by chance or depleted vs chance), then set to OR = 1 (pure chance)  
  OR[OR < 1] <- 1
  OR[pval >=sign.lev] <- 1
  
  #### Take log2 (e.g. if by chance OR = 1 => log2OR = 0. Also if log2OR = 1 => twice more common colocation
  #### then by chance)  
  OR <- log2(OR)
  #### Annotate with a STAR if positive OR is statistically significant e.g. p < 0.01)
  pval[OR < 1] <- 1
  pval[pval >=sign.lev] <- NA
  pval[pval < sign.lev] <- "*"
  cmap <- colorRampPalette(brewer.pal(9,"Blues"))(256)
  #Build heatmap based on log2OR correlation and annotate entry if OR is stat significant
  hm1 <- heatmap.2(OR,col=cmap,tracecol = NA,cellnote = pval,notecol = "black", main=txtstr1,key.title = "Color Key", key.xlab = "log2OR", key.ylab = "")
  if(outputval==1){
    return(hm1)}
  else{
    
    ###fisher output 
    
    cln<-c(1:(ncol(fishertable1)+1))
    pCol<-data.frame(t(c("pVal",rep("",2,ncol(fishertable1)))))
    odCol<-data.frame(t(c("Oddsratio",rep("",2,ncol(fishertable1)))))
    fishertable2<-data.frame(cbind(rownames(fishertable1),as.matrix(fishertable1)))
    fishertableodds2<-data.frame(cbind(rownames(fishertableodds1),as.matrix(fishertableodds1)))
    pName<-data.frame(t(colnames(fishertable2)))
    
    colnames(fishertableodds2)<-cln
    colnames(fishertable2)<-cln
    colnames(pName)<-cln
    colnames(pCol)<-cln
    colnames(odCol)<-cln
    
    out.df<-data.frame(rbind(pCol,pName,fishertable2,odCol,pName,fishertableodds2))
    return(out.df)
  }
}
heatmap_func<-function(datain,txtstr,outputval){
  ####Generates a normalized heatmap of all the genes and TFs in the selected GO-terms
  cmap <- colorRampPalette(brewer.pal(9,"Blues"))(256)
  data<-datain
  datamap<-as.matrix(data[,-1])
  datamap[is.infinite(datamap)]<-0
  datamap[datamap>1]<-1
  p1<-heatmap.2(datamap,col=cmap, tracecol=NA, main=txtstr)
  if(outputval==1){
    return(p1)}
  else{
    data2<-cbind(row.names(data),data[,-1])
    data2name<-c("Gene",colnames(data[,-1]))
    data2<-rbind(data2name,data2)
    out.df<-data2
    return(out.df)}
}
net_func<-function(datain, txtstr){
  ####Generates the network plots
  netMatrix <- as.matrix(datain)
  rownames(netMatrix) <- rownames(datain)
  netMatrix[netMatrix>1]<-1
  netMatrix[is.na(netMatrix)] <- 0
  dimnames(netMatrix) <- list(Gene = rownames(netMatrix),
                              TF = colnames(netMatrix))
  
  tfNetwork <- data.frame(netMatrix,
                          row.names = rownames(netMatrix))
  tfNetwork <- network(tfNetwork,
                       matrix.type = "bipartite",
                       ignore.eval = FALSE,
                       names.eval  = "weights")
  col <- c("Gene" = "#BAD0F0", "TF" = "#3F7DD3")
  tfNetwork %v% "Legend" = ifelse(network.vertex.names(tfNetwork) %in% rownames(netMatrix), "Gene", "TF")
  ggnet2(tfNetwork, alpha = 0.8, color.legend="Legend",size.legend="Connections",color.palette=col,label.size=3,label.color = "#404449", node.size="cent", node.color="Legend", size="degree",size.min=1,node.label = colnames(netMatrix), layout.par = list(niter = 1000), edge.size = 0.1)+ 
    guides(size = FALSE)+
    labs(title=txtstr)+theme_void()
  
}
model_zero_func<-function(datain, txtstr, outputval,name.TF){
  ####Generates the linear model for the selected GO-terms
  modTF="y~"
  for(i in 1:length(name.TF)){
    modTF<-cbind(modTF,paste(name.TF[i],"+", sep=""))}
  modTF<-paste(modTF,collapse="")
  modTF<-substr(modTF,1,nchar(modTF)-1)
  
  
  colnames(datain)[colnames(datain)=="Treated.Data$y"] <- "y"
  lm.test<-lm(modTF,datain)
  #dev.off()
  
  mlmrsum<-summary(lm.test)
  
  lm.test.fit<-data.frame(Fitted=fitted(lm.test), y.obs=datain$y)
  txtstr1<-capture.output(cat('r-sqr: ',mlmrsum$r.squared , 'r-sqradj: ',mlmrsum$adj.r.squared))
  
  p1<-ggplot(lm.test.fit, aes(Fitted, y.obs, group=Fitted)) +
    geom_smooth(method ="lm", color="blue", group=1)+
    geom_point()+
    geom_text(aes(label=rownames(datain)),hjust=0, vjust=0)+
    theme_classic()+
    labs(x="Fitted",y="log2(RPKM)", title=txtstr, subtitle=txtstr1)+
    theme(axis.text.y=element_text(size=10, colour="#4c4c4c", angle = 0, debug = FALSE), 
          axis.text.x = element_text(colour = "#4c4c4c", size=10, angle = 0, debug = FALSE),
          rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA))
  if( outputval==1){
    return(p1)}
  else{
    resmat<-data.frame(matrix(nrow=1,ncol=3))
    resmat<-cbind(data.frame(names(mlmrsum$residuals),(as.character(mlmrsum$residuals)),resmat))
    coefmat<-data.frame(mlmrsum$coefficients)
    coefmat<-cbind(rownames(coefmat),coefmat)
    coefmat<-rbind(c("","Estimate","Std. Error","t value","Pr(>|t|)"),coefmat)
    colnames(resmat)<-c(1:5)
    colnames(coefmat)<-c(1:5)
    coefname<-data.frame("Coefficients","","","","")
    colnames(coefname)<-c(1:5)
    resname<-data.frame("Residuals","","","","")
    colnames(resname)<-c(1:5)
    out.df<-data.frame(rbind(coefname,coefmat,resname,resmat))
    return(out.df)
  }
}
cluster_func<-function(datain,outputval,nclu){
  ####Use Partioning Around Medoids for the selected GO-terms using the number of peaks detected for each TF and each gene 
  data<-datain
  xclus<-as.matrix(data)
  xclus[which(!is.finite(xclus))] <- 0
  mydata <- na.omit(xclus) 
  mydata <- scale(xclus)
  mydata <- mydata[,colSums(is.na(mydata))<nrow(mydata)]
  pam.mydata= pam(mydata, nclu)
  p1<-fviz_cluster(pam.mydata, ggtheme=theme_pubr(),main=NULL)
  if(outputval==1){return(p1)}
  else{
    #Coef for the TFs in each cluster
    clusters<-pam.mydata$medoids
    clusters<-cbind(as.character(1:nclu),clusters)
    clusters<-rbind(colnames(clusters),clusters)
    clusters[1,1]<-"Cluster"
    colnames(clusters)<-c(1:(ncol(clusters)))
    #Which gene belongs to which cluster
    clustersdata<-as.data.frame(pam.mydata$clustering)
    clustersdata$names<-rownames(clustersdata)
    clustersdata<-rbind(c("Cluster", "Gene"),clustersdata)
    clusmat<-matrix(nrow=1,ncol=(ncol(clusters)-2))
    clustersdata<-cbind(clustersdata,clusmat)
    colnames(clustersdata)<-c(1:(ncol(clusters)))
    out.df<-data.frame(rbind(clusters,clustersdata))
    return(out.df)
  }
}


###Plotfunction for read profiles
alphaval=0.7
profile_plot_func<-function(TATA.df,temp_cond,x,Curr.Geneseq,ranges,colors.TF, data.BS, ATGC, input,yRanges,Gene_start,TPM, motiffinder){
 
  start.h<- max(3,temp_cond$y)
  if(input$yranges==TRUE){
    start.h<- max(yRanges)
  }
  p1 <- ggplot()+
    #Tatabox
  {if(input$TATA)geom_text(aes(x=TATA.df[1,1],y=1*start.h*0.03,label=TATA.df[1,3]),hjust=0, vjust=-0.5)}+
  {if(input$TATA)annotate("rect", xmin=TATA.df[1,]$x, xmax=TATA.df[1,]$y, ymin=(-1)*start.h*0.03 , ymax=1*start.h*0.03, alpha=alphaval, fill="Black" )}+
    #TF binding
    geom_line(data = temp_cond, aes(x = x, y = y, color = temp_cond$TF),size=1)+
    labs(colour = "TF")+
    ylab("norm. reads")+xlab("")+
    #Display Sequence
    {if(is.null(ranges$x[2])==TRUE)scale_x_continuous(breaks=x, labels=NULL)}+
    {if(sum(ranges$x[2]-ranges$x[1]>200))scale_x_continuous(breaks=x, labels=NULL)}+
    {if(sum(ranges$x[2]-ranges$x[1]<200))scale_x_continuous(breaks=x, labels=Curr.Geneseq)}+
    #Change the scales if user sppecify so
    coord_cartesian(xlim = ranges$x, expand = FALSE)+
    {if(input$yranges)coord_cartesian(ylim=yRanges,xlim = ranges$x, expand = FALSE)}+
    #CDS start
    geom_text(aes(x=Gene_start$x1,y=(start.h*0.03),label="CDS",hjust=0, vjust=-0.5))+
    annotate("rect", xmin=Gene_start$x1, xmax=Gene_start$x2, ymin=(-1)*start.h*0.03 , ymax=1*start.h*0.03, alpha=alphaval, fill="#5833c6" )+
    geom_text(aes(x=1900, y=start.h*0.8,label=paste("TPM",TPM) ))+
    #TSS
    geom_vline(xintercept = 1000,linetype="dotted")+
    geom_hline(yintercept = 0)+
    geom_text(aes(x=1002, y=start.h*0.7,label="TSS"), color="black",hjust=0)+
    #Set the colors legend for each TF and how many should be in ech column
    scale_color_manual(values=colors.TF)+
    guides(col = guide_legend(ncol = 3, nrow=8))+
    #Include GEM identified BS
    {if(input$TF_BS) annotate("rect", xmin=data.BS$x1, xmax=data.BS$x2, ymin=data.BS$y1*start.h*0.03 , ymax=data.BS$y2*start.h*0.03, alpha=alphaval, fill=data.BS$colval )}+
    #Allow searching for motifs
    {if(input$motiffinder)annotate("rect", xmin=motiffinder$x1, xmax=motiffinder$x2, ymin=motiffinder$y1*start.h*0.03 , ymax=motiffinder$y2*start.h*0.03, alpha=alphaval, fill=motiffinder$colval )}+
    theme(
      axis.text.x = element_text(colour = ATGC, hjust=1, angle = 0, debug = FALSE), #Sepcify the ATGC color code
      axis.text.y=element_text(size=10, colour="#4c4c4c", angle = 0, debug = FALSE),
      axis.line.x=element_line(color = NA),
      rect = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank())
}


function(input, output) {
  rdata.New.TF<-reactive({
        if(!input$submit1){
          for(cond in name.Cond){
            out_list<-c()
            WigLike<-NULL
            Wigname<-paste("wigLike",cond,sep="")
            out_list[Wigname]<-list(WigLike)
            out_list["GEManalysis"]<-NULL
            out_list["Targets"]<-NULL}
          
        }else{
          out_list<-c()
          withProgress(message = 'Uploading data',value = 0.1, {
          names<-c()
          name.file<-input$files[,1]
          for(cond in name.Cond){
          nameidx<-grep(paste("(.*)",cond,"(.*)wigLike$",sep=""),name.file)
          
          WigLike<-tryCatch({read.csv(input$files[[nameidx, 'datapath']],sep="\t",header=FALSE)},
                       error=function(file){
                         dummyvec<-data.frame(NA,NA,NA)
                         return(dummyvec) })
          
          names<-c(names,name.file[nameidx], sep="\n")
          Wigname<-paste("WigLike",cond,sep="")
          out_list[Wigname]<-list(WigLike)
          incProgress(0.7/length(name.Cond))
          }
          nameidx<-grep("(.*)_GEManalysis_(.*).csv$",name.file)
            
            names<-c(names,name.file[nameidx], sep="\n")
            GEManalysis <-tryCatch({read.csv(input$files[[nameidx, 'datapath']],sep=",")},
                                 error=function(file){
                                   dummyvec<-data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
                                   colnames(dummyvec)<-c("X","Gene","Chr","Pos","Strand","DistanceTSS","Strength","Condition","TF")
                                   return(dummyvec) })
          incProgress(0.1)
          nameidx<-grep("(.*)Target(.*)csv$",name.file)
          names<-c(names,name.file[nameidx], sep="\n")
          Targets <-tryCatch({read.csv(input$files[[nameidx, 'datapath']],sep=",")},
                                error=function(file){
                                  dummyvec<-data.frame(t(rep(NA, (length(name.Cond)+1))))
                                  colnames(dummyvec)<-c("X",name.Cond)
                                  return(dummyvec) })
          incProgress(0.1)
          output$text3<-renderText({names})
          out_list["GEManalysis"]<-list(GEManalysis)
          out_list["Targets"]<-list(Targets)
          })}
        return(out_list)
        })
  observeEvent(input$submit1,{
    call.data<-rdata.New.TF()
    output$text2<-renderText({"Data uploaded succesfully"})
  })
  rname.TF<-reactive({
    name.TF<-name.TF
    if(!input$submit1){name.TF<-name.TF}else{
      name.TF<-c(name.TF,"New_TF")
    }
    return(name.TF)
  })
  rmyColors.TF<-reactive({
    if(!input$submit1){ myColors.TF<-myColors.TF}else{
      myColors.TF<-c(myColors.TF,"black")
    }
    return(myColors.TF)
  })
  rdatafile_Peaks<-reactive({
    if(!input$submit1){
      datafile_Peaks<-datafile_Peaks}else{
        Gemanalysis<-rdata.New.TF()$GEManalysis
        Gemanalysis$TF<-"New_TF"
      datafile_Peaks<-rbind(datafile_Peaks,Gemanalysis)
      }
    return(datafile_Peaks)
    })

  output$Checkbox1<-renderUI({checkboxGroupInput("checkboxgrp1", "TFs", TF_list[1:(length(TF_list)/2)])})
  output$Checkbox2<-renderUI({checkboxGroupInput("checkboxgrp2", "TFs", TF_list[(length(TF_list)/2+1):length(TF_list)])})
  
  output$radioOut <- renderUI({
    if(input$submit1){
      new_list<-TF_list
      new_list[["New_TF"]]<-"New_TF"
      radioButtons("TFs","",new_list)}else{return(radioButtons("TFs","",TF_list))}
  })
  output$checkTF <- renderUI({
    if(input$submit1){
      checkboxInput("New_TF","New_TF",value=T)}else{return(NULL)}
  })
  output$text<-renderUI({
    geneS<-gene.Systematic
    geneC<-gene.Common
    colnames(geneS)<-"Systematic"
    colnames(geneC)<-"Common"
    autocomplete_list<-cbind(geneS,geneC)
    selectizeInput(inputId = "text",
                   label = "Select Gene",
                   choices = autocomplete_list,
                   selected = "GDH3",
                   multiple = FALSE, 
                   options = list(create = FALSE)) 
  })
  New_TF<- reactive({
    if(input$submit1){
      New_TF<-TRUE}else{New_TF<-FALSE}
  })
  rinputs<-reactive({
    activeTF<-c(input$checkboxgrp1[],input$checkboxgrp2[])
    inputs<-c()
    for(i in name.TF){
      inputs[i]<-FALSE
      if(!is.null(activeTF)){
        for(o in 1:length(activeTF)){
            if( as.character(activeTF[o])==as.character(i)){
              inputs[i]<-TRUE
            }
        }
      }
    }
     inputs["New_TF"]<-New_TF()
     return(inputs)
     })
  ranges <- reactiveValues(x = NULL, y = NULL)
  rGeneidx<-eventReactive(input$Load,{
    if(!input$text==""){
    input.text<-toupper(input$text)
    Geneidx<-grep(paste("^",input.text,"$",sep=""), as.character(gene.Common[,1]))
    if (length(Geneidx)==0){
      Geneidx<-grep(paste("^",input.text,"$",sep=""), as.character(gene.Systematic[,1]))
    }
    return(Geneidx)
    }else{return(NULL)}
  })

  rpeak.dist.data<-eventReactive(input$Load, {
   inputs<-rinputs()
   peak.dist.data<-c()
   for (cond in name.Cond){
      for (o in 1:length(rname.TF())){
        if(inputs[o]==TRUE){
          p  <-Reads_func(rname.TF()[o],cond,rdata.New.TF())
          lineDataCond=matrix(0L,nrow=1, ncol=2000)
          lineDataCond[1,p[p[,1]==as.character(gene.Common[rGeneidx(),1]),2]+1]<-p[p[,1]==as.character(gene.Common[rGeneidx(),1]),3]
          #Take every third data point
          y<-lineDataCond[,seq(1,2000,3)]
          data.dist<-data.frame(x.three,y,rname.TF()[o])
          colnames(data.dist)<-c("x","y","TF")
        }else{data.dist=(NULL)}
        if (o==1){
          temp.dist<-data.dist
        }else{
          temp.dist<-rbind(temp.dist, data.dist)
        }
      }

      peak.dist.data[cond]<-list(temp.dist)
    }
    return(peak.dist.data)
    })
  rBS<-eventReactive(input$Load,{
      BS<-c()
      input.gene<-as.character(gene.Systematic[rGeneidx(),1])
      input.gene2<-as.character(gene.Common[rGeneidx(),1])
      for(cond in name.Cond){
        datafile_Peaks<-rdatafile_Peaks()
        tempfile3<-subset(datafile_Peaks, datafile_Peaks$Gene == input.gene2 & datafile_Peaks$Condition == cond )
       
        x1<-c()
        x2<-c()
        y1<-c()
        y2<-c()
        t<-c()
        colval<-c()
        #Check if tempfile3 isn't 0
        if(nrow(tempfile3)>0){
          pos<-1000+as.numeric(tempfile3$DistanceTSS)
          #Generate the vector for plotting 
          for(j in 1:(length(pos))){
            x1[j]<-pos[j]-5
            x2[j]<-pos[j]+5
            y1[j]<-0.1
            y2[j]<-(-3.1)
            t[j]<-as.character(tempfile3$TF[j])
            colval[j]<-as.character(subset(rmyColors.TF(), rname.TF() == as.character(tempfile3$TF[j])))
          }
        }else{
          x1<-0
          x2<-0
          y1<-0
          y2<-0
          t<-"Seq"
          colval<-"#525252"}
        d=data.frame(x1, x2, y1, y2, t, colval)
        #Removes a BS if it is further than 1 or 2000 bp
        d<-subset(d, d$x1 > -1 & d$x2<2001)
        BS[cond]<-list(d)
      }
      return(BS)
  })
  rCurr.Geneseq<-reactive({
    if(is.null(ranges$x[1])==TRUE || sum(ranges$x[2]-ranges$x[1])>200){
      Curr.Geneseq=NULL}else{
        Geneidx<-rGeneidx()
        curr.geneSys<-gene.Systematic[Geneidx,1]
        Geneidx.seq<-grep(paste("^",curr.geneSys,"$",sep=""), as.character(Geneseq[,1]))
        Curr.Geneseq<-Geneseq[Geneidx.seq,-(1:2)]   
        Curr.Geneseq<-as.character(Curr.Geneseq)
        Curr.Geneseq<-sapply(seq(from=1, to=nchar(Curr.Geneseq), by=1), function(i) substr(Curr.Geneseq, i, i))
      }
    Curr.Geneseq
  })
  rATGC<-reactive({
    
    Curr.Geneseq<-rCurr.Geneseq()
    ATGC<-c()
    ATGC[Curr.Geneseq=="A"]<-("#f90000")
    ATGC[Curr.Geneseq=="T"]<-("#008e2a")
    ATGC[Curr.Geneseq=="G"]<-("#e2ad00")
    ATGC[Curr.Geneseq=="C"]<-("#113cfc")
    if(!is.null(ranges$x)){
      #xmin<-round((ranges$x[1]+0.1),0)
      #xmax<-round((ranges$x[2]+0.1),0)
      ATGC<-ATGC[ranges$x[1]:ranges$x[2]]
    }
    ATGC
  })
  ryRanges<-reactive({
    if(input$yranges==TRUE){
      max.temp<-c()
      for(cond in name.Cond){
        max.temp[cond]<-max(rpeak.dist.data()[[cond]]["y"])
      }
      max.yRanges<-max(max.temp[])
      if(is.infinite(max.yRanges) ==TRUE){
        yRanges <-c(-0.1,1)
        names(yRanges)<-c("ymin","ymax")
        return(yRanges )
      }else if (is.finite(max.yRanges)==TRUE){
        yRanges <-c(-max.yRanges*0.1,max.yRanges)
        names(yRanges)<-c("ymin","ymax")
        return(yRanges)
      }
    }
  })
  rTATA<-eventReactive(input$Load,{
      input.text<-toupper(input$text)
      Geneidx.tata<-grep(paste("^",input.text,"$",sep=""), as.character(TATA.file$gene_id) )
      if (length(Geneidx.tata)==0){
        Geneidx.tata<-grep(paste("^",input.text,"$",sep=""), as.character(TATA.file$gene_common))
      }
      
      TATAbox<-TATA.file[Geneidx.tata,"TATA_cen_pos"]
      if(is.na(TATAbox)==FALSE){
        TATAbox.x<-TATAbox-0.5
        TATAbox.y<-TATAbox+7.5
      }else if(is.na(TATAbox)==TRUE){
        TATAbox.x<-(1)*NA
        TATAbox.y<-(1)*NA
      }
      class<-c("TATAbox")
      colour<-c("grey")
      df<-data.frame(x=TATAbox.x, y=TATAbox.y,class=class,colour=colour)
      TATA_GTF.df<-df
  })
  
  rGene_start<-reactive({
    
    input.gene<-as.character(gene.Common[rGeneidx(),1])
    tempfile_start<-subset(TSS.start, TSS.start$Gene == input.gene)
    x1<-c(as.numeric(as.character(tempfile_start$DistanceTSStoORF))+1000)
    x2<-c(2000)
    y1<-c(-1)
    y2<-c(1)
    t<-c("start")
    colval<-c("#5833c6")
    d=data.frame(x1, x2, y1, y2, t, colval)
  })
  rTPM<-reactive({
    TPM<-c()
    input.gene<-as.character(gene.Systematic[rGeneidx(),1])
    tempfile.tpm<-subset(TPM.data, TPM.data$Gene == input.gene)
    
    for(cond in name.Cond){
     TPM[cond]<-data.frame(tempfile.tpm[cond])
     }
    TPM
  })
  rMotiffinder<-eventReactive(input$motiffinder,{
      Geneidx<-rGeneidx()
      curr.geneSys<-gene.Systematic[Geneidx,1]
      Geneidx.seq<-grep(paste("^",curr.geneSys,"$",sep=""), as.character(Geneseq[,1]))
      Curr.Geneseq<-Geneseq[Geneidx.seq,-(1:2)]   
      Curr.Geneseq<-as.character(Curr.Geneseq)
      motifstr<-toupper(input$motifstr)
      
      motifstr<-strsplit(motifstr, "[+]")
      motifstr<-unlist(motifstr)
      for(j in 1:length(motifstr)){
        
        dna.inputstr<-DNAString(motifstr[j])
        dna.Curr.Geneseq<-DNAString(Curr.Geneseq)
        motifposFW<-matchPattern(dna.inputstr,dna.Curr.Geneseq,max.mismatch=0,fixed="subject")
        dna.inputstr.rev<-reverseComplement(dna.inputstr)
        motifposRV<-matchPattern(dna.inputstr.rev,dna.Curr.Geneseq,max.mismatch=0,fixed="subject")
        pos.FW<-ranges(motifposFW)  
        pos.RV<-ranges(motifposRV)
        x1<-c()
        x2<-c()
        y1<-c()
        y2<-c()
        colval<-c()
        if(is.null(length(start(pos.FW)))&&is.null(length(start(pos.RV))))return(null)
        if(!is.null(length(start(pos.FW)))){
          for (i in 1:length(start(pos.FW))){
            x1[i]<-c(start(pos.FW)[i])
            x2[i]<-c(end(pos.FW)[i])
            y1[i]<-c(1)
            y2[i]<-c(-1)
            colval[i]<-c("blue")
          }
          d.FW<-data.frame(x1, x2, y1, y2, colval)
        }
        x1<-c()
        x2<-c()
        y1<-c()
        y2<-c()
        colval<-c()
        if(!is.null(length(start(pos.RV)))){
          for (i in 1:length(start(pos.RV))){
            x1[i]<-c(start(pos.RV)[i])
            x2[i]<-c(end(pos.RV)[i])
            y1[i]<-c(1)
            y2[i]<-c(-1)
            colval[i]<-c("red")
          }
          d.RV<-data.frame(x1, x2, y1, y2, colval)}
        if (j==1){
          d<-rbind(d.FW,d.RV)  
        }
        if(j>1){
          d<-rbind(d,d.FW,d.RV)
        }
        return(d)}
  })
  rdataTarget<-reactive({

  if(!is.null(input$TFs)){

    if(input$TFs=="New_TF"){ 
        dataTarget<-rdata.New.TF()$Targets
    }else{
      
      name<-list.files(path="TF_data_files/Data/",pattern=paste(input$TFs,"_geneTargetList_(.*).csv$",sep=""))
      dataTarget <-tryCatch({read.csv(paste0("TF_data_files/Data/",name),sep=",")},
                          error=function(file){
                            dummyvec<-data.frame(t(rep(NA, (length(name.Cond)+1))))
                            colnames(dummyvec)<-c("X",name.Cond)
                            return(dummyvec) })
      }
      colnames(dataTarget)[1]<-"Gene Common"
         DT::datatable(dataTarget, rownames=FALSE, extensions = c('FixedColumns',"FixedHeader"), 
                  options = list(dom = 't', 
                                 scrollX = TRUE, 
                                 paging=FALSE,
                                 fixedHeader=TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
    
      }else{dataTarget=NULL}
    return(dataTarget)
  })
  
  observeEvent(input$plot1_dblclick, {
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(round(brush$xmin+0.1,0), round(brush$xmax+0.1,0))
      names(ranges$x)<-c("xmin","xmax")
      ranges$y <- c(brush$ymin, brush$ymax)
      names(ranges$y)<-c("ymin","ymax")
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 0),  "\n")
    }
    
    xcord<-round(as.numeric(input$plot1_hover[1]),0)
    paste0(
      "Distance to TSS: ", (xcord-1000)
    )})
  output$GeneInfo<-renderText({
    if(length(rGeneidx())>0){
      curr.geneSys<-gene.Systematic[rGeneidx(),1]
      input.gene<-as.character(gene.Common[rGeneidx(),1])
      tempfile_start<-subset(TSS.start, TSS.start$Gene == input.gene)
      Strand<-tempfile_start$Strand
      Chr<-tempfile_start$Chromosome
      ChrPos<-as.numeric(as.character(tempfile_start$TSS))+as.numeric(as.character(tempfile_start$DistanceTSStoORF))
      if(input.gene!="GDH3" & input.gene!="ARR1"){
        nextgeneS<-gene.Systematic[(rGeneidx()+1),1]
        prevgeneS<-gene.Systematic[(rGeneidx()-1),1]
        nextgeneC<-gene.Common[(rGeneidx()+1),1]
        prevgeneC<-gene.Common[(rGeneidx()-1),1]
        if(substr(curr.geneSys,2,2)!=substr(nextgeneS,2,2)){
          nextgeneC<-"End of Chromosome"
          nextgeneS<-""}
        if(substr(curr.geneSys,2,2)!=substr(prevgeneS,2,2)){
          prevgeneC<-"End of Chromosome"
          prevgeneS<-""
        }}else if(input.gene=="GDH3" ){
          nextgeneS<-gene.Systematic[(rGeneidx()+1),1]
          nextgeneC<-gene.Common[(rGeneidx()+1),1]
          prevgeneC<-"End of Chromosome"
          prevgeneS<-"" 
        }else if(input.gene=="ARR1" ){
          prevgeneS<-gene.Systematic[(rGeneidx()-1),1]
          prevgeneC<-gene.Common[(rGeneidx()-1),1]
          nextgeneC<-"End of Chromosome"
          nextgeneS<-""
        }
      
      return(paste(paste("Genename Systematic: ", curr.geneSys), paste("Genename Common: ", input.gene),paste("Gene start: ",Chr,ChrPos),paste("Gene orientation: ",Strand),paste("Gene up: ", prevgeneC,prevgeneS),paste("Gene down: ", nextgeneC, nextgeneS),sep="\n"))
    }
    else{return("Error: Gene doesn't exist")}
  })
  rSeqMax<-eventReactive(input$SeqFind,{
    SeqMax<-input$SeqMax
  })
  rSeqMin<-eventReactive(input$SeqFind,{
    SeqMin<-input$SeqMin
  })
  
  output$Sequnce_out <- renderText({
      SeqMax<-rSeqMax()
      SeqMin<-rSeqMin()
    if((as.numeric(SeqMax)-as.numeric(SeqMin))>0){
      curr.geneSys<-gene.Systematic[rGeneidx(),1]
      Geneidx.seq<-grep(paste("^",curr.geneSys,"$",sep=""), as.character(Geneseq[,1]))
      Curr.Geneseq<-as.character(Geneseq[Geneidx.seq,-(1:2)] )  
      Curr.Geneseq<-sapply(seq(from=1, to=nchar(Curr.Geneseq), by=1), function(i) substr(Curr.Geneseq, i, i))
      Curr.Geneseq<-Curr.Geneseq[1000+(SeqMin:SeqMax)]
      d<-floor(((1000+as.numeric(SeqMax))-(1000+as.numeric(SeqMin)))/20)
      Out.Geneseq<-data.frame(matrix(ncol=d,nrow=1))
      o=0
      if((as.numeric(SeqMax)-as.numeric(SeqMin))<20){
        Out.Geneseq<-Curr.Geneseq
      }else{
        for (i in 1:d){
          if(i>1){
            o=1}
          Out.Geneseq[,i]<-paste(Curr.Geneseq[((i-1)*20+o):(i*20)],collapse="")
        }
      }
      Out.Geneseq<-as.character(Out.Geneseq)
      
    }else(Out.Geneseq<-"Wrong input sequence")
      return(Out.Geneseq)
    })

  rCCM<-reactive({
    val1=0
    if(input$val1=="Exclude dubious and non-Yeast 8" ){val1=1}
    if(input$val1=="Include dubious and non-Yeast 8"){val1=2}
    return(val1)
  })
  rGOterm<-eventReactive(input$search,{
    goterm<-input$goterm
  })
  rTest<-eventReactive(input$search,{
    test<-input$test
  })
  
  output$goterms<-renderDataTable({
    pathways.name<-unlist(strsplit(tolower(input$goterm), "[+]"))
    pathways.go <- filter(go_bio, grepl(paste(pathways.name, collapse="|"),name_1006))
    pathways.go1<-data.frame(unique(pathways.go$name_1006))
    colnames(pathways.go1)<-"GO-terms"
    return(pathways.go1)
  }, options=list(pageLength=5))
  
  output$TestInput<-renderUI({
    if(input$test=="Fisher"){
      textOut<-("The Fisher test is based on the occurrence of one transcription factor at the same gene as another transcription factor. The probability is plotted in a heatmap where * denotes p<0.001")
    }
    if(input$test=="Heatmap"){
      textOut<-("Generates a heatmap of the selected genes where the total number of binding sites are normalized to 1")
    }
    if(input$test=="Network"){
      textOut<-("The Network shows each transcription factor and its targets vs all other transcription factors. Transcription factors that are highly connected appear closer on the chart. The node size is weighted based on the number of edges.")
    }
    if(input$test=="Cluster"){
      textOut<-("The genes in the GO-term is clustered by k-medoids clustering based on the transcription factor occurrence. Positive medoid coefficient indicates presence of TF while negative coefficient indicates absence of TF in the cluster.")
    }
    if(input$test=="Linear Model"){
      textOut<-("Based on the transcription factor occurrence a linear model is used to predict the resulting TPM values. Both the R-square and the R-squared adjusted is displayed" )
    }
    tagList(
      actionButton("TestInfo", "?"),
      bsTooltip("TestInfo", textOut,
                "right", options = list(container = "body"))
    )
  })

  output$table1 <- renderDataTable({ rdataTarget()}, options=list(pageLength=5))
  output$sctrex<-renderImage({
    
    return(list(
      src = "TF_data_files/Resources/sctrex.png",
      contentType = "image/png",
      alt = "sctrex",
      height    = 200,
      units     = "in"
    ))
    
  },deleteFile = FALSE)

  output$plots.peak <- renderUI({
    plot_output_list <- lapply(1:length(name.Cond), function(i) {
      plotname <- paste("plot", i, sep="")
      plottitle <- paste("plottitle", i, sep="")
      tags$div(class = "group-output",
               textOutput(plottitle, container = h3),
               plotOutput(plotname, 
                          height = 180, 
                          dblclick = "plot1_dblclick", 
                           hover ="plot1_hover",
                           brush = brushOpts(id = "plot1_brush",resetOnNew = TRUE)
               ) %>% withSpinner(color="#0dc5c1")
               
      )
    })
    do.call(tagList, plot_output_list)
  })
  output$plots.analysis <- renderUI({
    analysis_output_list <- lapply(1:length(name.Cond), function(i) {
      plotname <- paste("plotAn", i, sep="")
      plottitle <- paste("plottitleAn", i, sep="")
      tags$div(class = "group-output",
               textOutput(plottitle, container = h3),
               plotOutput(plotname, 
                          height = 350 
               ) %>% withSpinner(color="#0dc5c1")
               
      )
    })
    do.call(tagList, analysis_output_list)
  })
  output$plots.motif <- renderUI({
    motif_output_list <- lapply(1:length(name.Cond), function(i) {
      plotname <- paste("Motif", i, sep="")
      plottitle <- paste("plottitleMot", i, sep="")
      tags$div(class = "group-output",
               textOutput(plottitle, container = h3),
               imageOutput(plotname,height=200))
    })
    do.call(tagList, motif_output_list)
  })
  output$plots.seq <- renderUI({
    seq_output_list <- lapply(1:length(name.Cond), function(i) {
      plotname <- paste("Seq", i, sep="")
      plottitle <- paste("plottitleSeq", i, sep="")
      tags$div(class = "group-output",
               textOutput(plottitle, container = h3),
               imageOutput(plotname,height=200))
    })
    do.call(tagList, seq_output_list)
  })
  output$plots.PeakDist <- renderUI({
    PeakDist_output_list <- lapply(1:length(name.Cond), function(i) {
      plotname <- paste("PeakDist", i, sep="")
      plottitle <- paste("plottitlePeakDist", i, sep="")
      tags$div(class = "group-output",
               textOutput(plottitle, container = h3),
               imageOutput(plotname,height=200))
    })
    do.call(tagList, PeakDist_output_list)
  })
  output$plot5 <-renderPlot({
    n=length(rname.TF())
    
    n2=10
    x1<-matrix(0,1,n)
    x2<-matrix(0,1,n)
    y1<-matrix(0,1,n)
    y2<-matrix(0,1,n)
    t<-rname.TF()
    mycol<-as.character(rmyColors.TF())
    names(mycol)<-as.character(t(rname.TF()))
    colval<-mycol
    d=data.frame(t(x1), t(x2), t(y1), t(y2), t, colval)
    colnames(d)<-c("x1","x2","y1","y2","t","colval")
    p1<-ggplot() +
      scale_x_continuous(name="x") +
      scale_y_continuous(name="y")+
      labs(fill = "TF BS", alpha=0.5)+
      scale_fill_manual(values=mycol)+ 
      geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), alpha=alphaval) +
      theme_void()+
      guides(fill = guide_legend(nrow = nrow(rname.TF()), ncol=1))+
      theme(legend.text=element_text(size=14),
            legend.title=element_text(size=16),
            panel.grid.major = element_blank() # get rid of major grid
            , panel.grid.minor = element_blank() # get rid of minor grid
      )
    
    p1
  },bg="transparent")
  #Dynamically generate all the plots based on conditions
  for (i in 1:length(name.Cond)) {
    local({
      my_i<-i
      plottitleMot <- paste("plottitleMot", my_i, sep="")
      output[[plottitleMot]] <- renderText({name.Cond[my_i]})
      plotnameMot <- paste("Motif", my_i, sep="")
      output[[plotnameMot]] <- renderImage({
        motif<-motif_func(input$TFs,name.Cond[my_i])
      },deleteFile = FALSE)
      
      plottitleSeq <- paste("plottitleSeq", my_i, sep="")
      output[[plottitleSeq]] <- renderText({name.Cond[my_i]})
      plotnameSeq <- paste("Seq", my_i, sep="")
      output[[plotnameSeq]] <- renderImage({
        seq_map<-seq_map_func(input$TFs,name.Cond[my_i])
      },deleteFile = FALSE)
      
      plottitlePeakDist <- paste("plottitlePeakDist", my_i, sep="")
      output[[plottitlePeakDist]] <- renderText({name.Cond[my_i]})
      plotnamePeakDist <- paste("PeakDist", my_i, sep="")
      output[[plotnamePeakDist]] <- renderImage({
        peak_dist<-peak_dist_func(input$TFs,name.Cond[my_i])
      },deleteFile = FALSE)
      
      plottitle <- paste("plottitle", my_i, sep="")
      output[[plottitle]] <- renderText({name.Cond[my_i]})
      plotname <- paste("plot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        mycol<-as.character(rmyColors.TF())
        names(mycol)<-as.character(t(rname.TF()))
        if (length(rpeak.dist.data()[[name.Cond[my_i]]])==0)return(NULL)
        p1<-profile_plot_func(rTATA(),rpeak.dist.data()[[name.Cond[my_i]]], x,rCurr.Geneseq(),ranges,mycol, rBS()[[name.Cond[my_i]]], rATGC(), input, ryRanges(),rGene_start(),rTPM()[[name.Cond[my_i]]],rMotiffinder())
        p1

      },bg="transparent")
      
      plottitleAn <- paste("plottitleAn", my_i, sep="")
      output[[plottitleAn]] <- renderText({name.Cond[my_i]})
      plotnameAn <- paste("plotAn", my_i, sep="")
      output[[plotnameAn]] <- renderPlot({
        downloadVal=1
        p1<-tryCatch({Stat_plot_func(rCCM(),downloadVal,name.Cond[my_i],rGOterm(),rname.TF(),rTest(),rdata.New.TF(),input$slider1)},
                               error=function(analysis){
                                 return(ggplot()+theme_void()+geom_text(aes(x=1,y=1,label="Not enough genes selected to do statistical function",hjust=(0.4) ),size=5))}
        )
        p1
      },bg="transparent")
      
    })
  }
  
  output$ReadsDown<-renderUI({
    selectInput("ReadsData", "Download Data",name.Cond, selected=name.Cond[1], selectize=TRUE)
  })
  output$StatDown<-renderUI({
    selectInput("StatData", "Download Data",name.Cond, selected=name.Cond[1], selectize=TRUE)
  })
  output$downloadreadsData <- downloadHandler(
    filename = function() {
      paste(input$text,"_",input$ReadsData, "_Reads.csv", sep = "")
    },
    content = function(file) {
      
      peak.dist.data<-rpeak.dist.data()
      out<-peak.dist.data[input$ReadsData]
      write.csv(out, file, row.names = FALSE)
    }
  )
  output$downloadPeakData <- downloadHandler(
    filename = function() {
      paste(input$TFs, "_SumPeaks.csv", sep = "")
    },
    content = function(file) {
      file.name<-list.files(path="TF_data_files//",pattern=paste(input$TFs,"(.*)_geneTargetList_(.*).csv$",sep=""))
      out<-read.csv(paste("TF_data_files/Data/",file.name,sep=""),sep=",")
      write.csv(out, file, row.names = FALSE)
    }
  )
  output$downloadPeakDataAd <- downloadHandler(
    filename = function() {
      paste(input$TFs, "_All_PeakPos.csv", sep = "")
    },
    content = function(file) {
      file.name<-list.files(path="TF_data_files/Data/",pattern=paste(input$TFs,"(.*)_GEManalysis_(.*).csv$",sep=""))
      out<-read.csv(paste("TF_data_files/Data/",file.name,sep=""),sep=",")
    
      write.csv(out, file, row.names = FALSE)
    }
  )
  rdowloadStatData<-reactive({
    downloadVal=2
    out.data<-Stat_plot_func(rCCM(),downloadVal,input$StatData,rGOterm(),rname.TF(),rTest(),rdata.New.TF(),input$slider1)
    return(out.data)
  })
  output$downloadStatData <- downloadHandler(
    filename = function() {
      paste(input$StatData,"_",input$test, ".csv", sep = "")
    },
    content = function(file) {
      write.table(rdowloadStatData(), file, row.names = FALSE, col.names=FALSE, sep=",")
    }
  ) 
  
  }

