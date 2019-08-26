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
colorset.TF.df<-data.frame(name.TF,myColors.TF)



#Load TSS start position
TSS.start<-read.csv(paste("TF_data_files/Resources/","190704_TSSData.tsv",sep=""), sep="\t", header=TRUE)
#Get the current genelist
geneList<-data.frame(TSS.start$GeneName, TSS.start$Gene)
colnames(geneList)<-c("Gene","GeneC")
gene.Systematic<-data.frame(geneList[,1])
colnames(gene.Systematic)<-c("GeneS")
gene.Common<-data.frame(geneList[,2])
colnames(gene.Common)<-c("GeneC")
#Load metabolic genes from Yeast7
metabolicgenes<-read.csv(paste("TF_data_files/Resources/","Yeast8_genes.csv",sep=""),sep=";", header=TRUE)
#Load TATA and GTF positions. GTF is not used in this version 
TATA.GTF.file<-read.csv(paste("TF_data_files/Resources/","TATA_pos.csv",sep=""),sep=";", header=TRUE)
#Load Sequence for each Gene
Geneseq<-read.csv(paste("TF_data_files/Resources/","TSS_anotation_seq.bed",sep=""),sep="\t",header = FALSE )
#Load TPM
TPM_data<-read.csv(paste("TF_data_files/Resources/","TPM.csv",sep=""),sep=";")
#create dataframes for TFs
x  <- matrix(1:2000)
x.three<-seq(1,2000,3)

Reads_func<- function(tf,cond) {
  p = read.csv(paste("TF_data_files/WigLike/",tf,"_",cond,"_ol_combRep_geneAssigned_190314.wigLike",sep=""),header = FALSE,sep = "\t")
  }
#Load GEM peaks 
for (i in name.TF){
  if(i=="Cat8"){
    datafile<-read.csv(paste("TF_data_files/GEMPeaks/",i,"_GEManalysis_190314.csv",sep=""),sep=",")
    datafile$TF<-i
  }else{
    datafile2<-read.csv(paste("TF_data_files/GEMPeaks/",i,"_GEManalysis_190314.csv",sep=""),sep=",")
    datafile2$TF<-i
    datafile<-rbind(datafile,datafile2)}
  datafile_Peaks<-datafile
}

for (i in name.TF){
  if(i=="Cat8"){
    datafile<-read.csv(paste("TF_data_files/GEMPeaks/",i,"_geneTargetList_190314.csv",sep=""),sep=",")
    datafile$TF<-i
  }
  datafile2<-read.csv(paste("TF_data_files/GEMPeaks/",i,"_geneTargetList_190314.csv",sep=""),sep=",")
  datafile2$TF<-i
  datafile_SumPeaks<-rbind(datafile,datafile2)
}
#Generate data for statistical analysis 
for (c in name.Cond){
  data_out<-geneList
  for (i in name.TF){
    datafile<-read.csv(paste("TF_data_files/GEMPeaks/",i,"_geneTargetList_190314.csv",sep=""),sep=",")
    condTF<-datafile[names(datafile) == c] 
    datafile.2<-data.frame(datafile$X,condTF)
    colnames(datafile.2)<-c("GeneC",i)
    data_out<-merge(data_out,datafile.2, by="GeneC", all=TRUE)
  }
  assign(paste(c,"_Peaks", sep=""),data_out)
}




####Load GO-terms located in the gogenes file, if not then load from Ensembl
# mtry <- try(read.csv(paste("TF_data_files/Resources/","gogenes.csv",sep=""),sep=";"), 
#             silent = TRUE)
# if (class(mtry) != "try-error") {
#   go_bio<-read.csv(paste("TF_data_files/Resources/","gogenes.csv",sep=""),sep=";")
# } else {
  # message("File doesn't exist, loading from Ensembl")
  ensembl=useMart("ENSEMBL_MART_ENSEMBL", "scerevisiae_gene_ensembl", host="www.ensembl.org")
  gogenes<-getBM(attributes=c('ensembl_gene_id','name_1006','namespace_1003'),mart=ensembl)
  
  go_bio<-gogenes[grep("biological_process",gogenes$namespace_1003),] # Select only Biological Process
  go_bio<-suppressMessages(remove.vars(go_bio,c("namespace_1003")))
  
# }
go_bio_unique<-data.frame(unique(go_bio$name_1006))



###Functions for generating Fisherplots, Net plots, models or Clusters
goterm_func <- function(datain) {
  ####Locates the GO-terms included in the users searchterm
  pathways.name<-datain
  pathways.go <- filter(go_bio, grepl(paste(pathways.name, collapse="|"),name_1006))
  
  pathways.gene<-data.frame(unique(pathways.go$ensembl_gene_id))
  colnames(pathways.gene)<-"Gene"
  return(pathways.gene)
}
data_func<- function(datain,pathway.val,val1,val2) {
  ####Generates the data used, if Yeast8 and dubious is or isn't selected
  dataset<-c("_Peaks","_Peaks_SN","_Sum")
  data1<-get(paste(datain,dataset[val2], sep=""))
  ####Remove non-metabolic genes
  if (val1==1){
    data1<-merge(data1,metabolicgenes,by="Gene")
  }
  ####Merge with pathway genes
  data1<-merge(data1,pathway.val,by="Gene")
  if(is.null(data1)){data1<-c("Gene"="empty")}
  return(list("data1"=data1))
}
data_treatment_func <- function(datain1,datain2,datain3, val2, name.Cond) {
  y1<-datain1
  x1<-datain2
  rownames(x1)<-datain3$CommonName
  ####Only include genes which has a TPM above 1
  x1<-x1[y1>1,]
  x1[is.na(x1)]<-0
  y1<-y1[y1>1]
  y1<-log2(y1)
  if (val2==3){
    data1<-get(paste(name.Cond,"_Sum", sep=""))
    data1<-data1[complete.cases(data1), ]
    av<-mean(as.matrix(data1[,2:length(name.TF)]))
    x1[x1<(0.5*av)]<-0
  }
  ####only include genes which has at least one binding or binding higher than 50% average 
  y1<-y1[rowSums(x1)>0]
  x1<-x1[rowSums(x1)>0,]
  x1[is.na(x1)]<-0
  return(list("y"=y1,"x"=x1))
}
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
    fishertable1<-data.frame(cbind(rownames(fishertable1),fishertable1))
    fishertableodds1<-data.frame(cbind(rownames(fishertableodds1),as.matrix(fishertableodds1)))
    pName<-data.frame(t(colnames(fishertable1)))
    
    colnames(fishertableodds1)<-cln
    colnames(fishertable1)<-cln
    colnames(pName)<-cln
    colnames(pCol)<-cln
    colnames(odCol)<-cln
    
    out.df<-data.frame(rbind(pCol,pName,fishertable1,odCol,pName,fishertableodds1))
    return(out.df)
  }
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
  nodsiz = c("Gene" = 1, "TF" = 11)
  tfNetwork %v% "Legend" = ifelse(network.vertex.names(tfNetwork) %in% rownames(netMatrix), "Gene", "TF")
  ggnet2(tfNetwork, alpha = 0.8, color.legend="Legend",size.legend="Connections",color.palette=col,label.size=3,label.color = "#404449", node.size="cent", node.color="Legend", size="degree",size.min=1,node.label = colnames(netMatrix), layout.par = list(niter = 1000), edge.size = 0.1)+ 
    guides(size = FALSE)+
    labs(title=txtstr)+theme_void()
  
}
model_zero_func<-function(datain, txtstr, outputval){
  ####Generates the linear model for the selected GO-terms
  modTF="y~"
  for(i in 1:length(name.TF)){
    modTF<-cbind(modTF,paste(name.TF[i],"+", sep=""))}
  modTF<-paste(modTF,collapse="")
  modTF<-substr(modTF,1,nchar(modTF)-1)
  
  
  colnames(datain)[colnames(datain)=="data.list2$y"] <- "y"
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
  mydata <- na.omit(xclus) # listwise deletion of missing
  mydata <- scale(xclus) # standardize variables
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
alphaval=0.7

profile_plot_func<-function(TATA_GTF.df,temp_cond, x,Curr.Geneseq,ranges,colorset.TF, data.BS, ATGC, input,yRanges,Gene_start,TPM, motiffinder){
  
  start.h<- max(3,temp_cond$y)
  if(input$yranges==TRUE){
    start.h<- max(yRanges)
  }
  p1 <- ggplot()+
    #Tatabox
  {if(input$GTF)geom_text(aes(x=TATA_GTF.df[1,1],y=1*start.h*0.03,label=TATA_GTF.df[1,3]),hjust=0, vjust=-0.5)}+
  {if(input$GTF)annotate("rect", xmin=TATA_GTF.df[1,]$x, xmax=TATA_GTF.df[1,]$y, ymin=(-1)*start.h*0.03 , ymax=1*start.h*0.03, alpha=alphaval, fill="Black" )}+
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
    scale_color_manual(values=colorset.TF)+
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
  rinputs<-reactive({
    #Add TF
    inputs<-data.frame(input$Cat8, input$Cbf1,input$Ert1, input$Gcn4, input$Gcr1, input$Gcr2,input$Hap1,input$Ino2,input$Ino4,input$Leu3,input$Oaf1,input$Pip2,input$Rds2,input$Rgt1,input$Rtg1,input$Rtg3,input$Sip4,input$Stb5,input$Sut1,input$Tye7)
    inputs})
  ranges <- reactiveValues(x = NULL, y = NULL)
  rGeneidx<-eventReactive(input$Load,{
    input.text<-toupper(input$text)
    Geneidx<-grep(paste("^",input.text,"$",sep=""), as.character(gene.Common[,1]))
    if (length(Geneidx)==0){
      Geneidx<-grep(paste("^",input.text,"$",sep=""), as.character(gene.Systematic[,1]))
    }
    Geneidx
  })

  rpeak.dist.data<-eventReactive(input$Load, {
   inputs<-rinputs()
    for (i in 1:4){
      for (o in 1:length(name.TF)){
        if(inputs[o]==TRUE){
          p  <-Reads_func(name.TF[o],name.Cond[i])

          lineDataCond=matrix(0L,nrow=1, ncol=2000)
          lineDataCond[1,p[p[,1]==as.character(gene.Common[rGeneidx(),1]),2]+1]<-p[p[,1]==as.character(gene.Common[rGeneidx(),1]),3]
          #Take every third data point
          y<-lineDataCond[,seq(1,2000,3)]
          data.dist<-data.frame(x.three,y,name.TF[o])
          colnames(data.dist)<-c("x","y","TF")
        }else{data.dist=(NULL)}
        if (o==1){
          temp.dist<-data.dist
        }else{
          temp.dist<-rbind(temp.dist, data.dist)
        }
        
      }
      if(i==1){
        
        tflist_Glu<-temp.dist}
      if(i==2){
        
        tflist_Nit<-temp.dist}
      if(i==3){
        
        tflist_Eth<-temp.dist}
      if(i==4){
        
        tflist_Ana<-temp.dist}
    }
    return(peak.dist.data<-list("tflist_Glu"=tflist_Glu,"tflist_Eth"=tflist_Eth,"tflist_Nit"=tflist_Nit,"tflist_Ana"=tflist_Ana))
    })
  rBS<-reactive({
    
    if(input$TF_BS==TRUE){
      input.gene<-as.character(gene.Systematic[rGeneidx(),1])
      input.gene2<-as.character(gene.Common[rGeneidx(),1])
      for(i in 1:4){
        tempfile3<-subset(datafile_Peaks, datafile_Peaks$Gene == input.gene2 & datafile_Peaks$Condition == name.Cond[i] )#& datafile$S>4)
        x1<-c()
        x2<-c()
        y1<-c()
        y2<-c()
        t<-c()
        colval<-c()
        #Check if tempfile3 isn't 0
        if(nrow(tempfile3)>0){
          pos<-1000-as.numeric(tempfile3$DistanceTSS)
          #Generate the vector for plotting 
          for(j in 1:(length(pos))){
            x1[j]<-pos[j]-5
            x2[j]<-pos[j]+5
            y1[j]<-0.1
            y2[j]<-(-3.1)
            t[j]<-as.character(tempfile3$TF[j])
            colval[j]<-as.character(subset(colorset.TF.df$myColors, colorset.TF.df$name.TF == as.character(tempfile3$TF[j])))
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
        if(i==1){d.Glu<-d}
        if(i==2){d.Nit<-d}
        if(i==3){d.Eth<-d}
        if(i==4){d.Ana<-d}
      }
      d<-list("d.Glu"=d.Glu,"d.Eth"=d.Eth,"d.Nit"=d.Nit,"d.Ana"=d.Ana)
    }
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
      peak.dist.data<-rpeak.dist.data()
      max.temp.glu<-max(peak.dist.data$tflist_Glu["y"])
      max.temp.eth<-max(peak.dist.data$tflist_Eth["y"])
      max.temp.ana<-max(peak.dist.data$tflist_Ana["y"])
      max.temp.nit<-max(peak.dist.data$tflist_Nit["y"])
      
      max.yRanges<-max(c(max.temp.glu,max.temp.eth,max.temp.ana,max.temp.nit))
      
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
  rTATA_GTF<-reactive({
    if(input$GTF==TRUE){
      input.text<-toupper(input$text)
      # input.text<-"ENO1"
      Geneidx.tata<-grep(paste("^",input.text,"$",sep=""), as.character(TATA.GTF.file$gene_id) )
      if (length(Geneidx.tata)==0){
        Geneidx.tata<-grep(paste("^",input.text,"$",sep=""), as.character(TATA.GTF.file$gene_common))
      }
      
      TATAbox<-TATA.GTF.file[Geneidx.tata,"TATA_cen_pos"]
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
    }
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
    input.gene<-as.character(gene.Systematic[rGeneidx(),1])
    tempfile_tpm<-subset(TPM_data, TPM_data$Gene == input.gene)
    TPM<-data.frame(tempfile_tpm$Glu.lim,tempfile_tpm$Nit.lim,tempfile_tpm$Eth.lim,tempfile_tpm$Ana.lim)
    colnames(TPM)<-c("Glu","Nit","Eth","Ana")
    TPM
  })
  rMotiffinder<-reactive({
    if(input$motiffinder==TRUE){
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
      }
      
      return(d)
    }
  })
  rdataTarget<-reactive({
    dataTarget<-read.csv(paste("TF_data_files/GEMPeaks/",input$TFs,"_geneTargetList_190314.csv",sep=""),sep=",")
    colnames(dataTarget)<-c("Gene Common","Glu-lim","Nit-lim","Eth-lim","Ana-lim")
    DT::datatable(dataTarget, rownames=FALSE, extensions = c('FixedColumns',"FixedHeader"), 
                  options = list(dom = 't', 
                                 scrollX = TRUE, 
                                 paging=FALSE,
                                 fixedHeader=TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
    return(dataTarget)
  })
  
  output$table1 <- renderDataTable({ rdataTarget()}, options=list(pageLength=5))
  output$motif_glu<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Glu_MEME/logo1.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$motif_nit<-renderImage({
    TFconsensusname<-paste("TF_data_files//Logos/",input$TFs,"_Nit_MEME/logo1.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$motif_ana<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Ana_MEME/logo1.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$motif_eth<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Eth_MEME/logo1.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  
  output$motif_glu_seq<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Glu_MEME/",input$TFs,"_Glu_SeqMap.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$motif_nit_seq<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Nit_MEME/",input$TFs,"_Nit_SeqMap.png",sep="")

    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$motif_ana_seq<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Ana_MEME/",input$TFs,"_Ana_SeqMap.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    },deleteFile = FALSE)
  output$motif_eth_seq<-renderImage({
    TFconsensusname<-paste("TF_data_files/Logos/",input$TFs,"_Eth_MEME/",input$TFs,"_Eth_SeqMap.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "motif",
      height    = 100,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  
  output$peakdist_glu<-renderImage({
    TFconsensusname<-paste("TF_data_files/PeakDist/",input$TFs,"_Glu_PeakHistogram_190314.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "peakdist",
      height    = 170,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$peakdist_nit<-renderImage({
    TFconsensusname<-paste("TF_data_files/PeakDist/",input$TFs,"_Nit_PeakHistogram_190314.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "peakdist",
      height    = 170,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$peakdist_eth<-renderImage({
    TFconsensusname<-paste("TF_data_files/PeakDist/",input$TFs,"_Eth_PeakHistogram_190314.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "peakdist",
      height    = 170,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  output$peakdist_ana<-renderImage({
    TFconsensusname<-paste("TF_data_files/PeakDist/",input$TFs,"_Ana_PeakHistogram_190314.png",sep="")
    return(list(
      src = TFconsensusname,
      contentType = "image/png",
      alt = "peakdist",
      height    = 170,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  
  output$sctrex<-renderImage({
    
    return(list(
      src = "TF_data_files/Resources/sctrex.png",
      contentType = "image/png",
      alt = "sctrex",
      height    = 200,
      units     = "in"
    ))
    
  },deleteFile = FALSE)
  
  output$downloadreadsData <- downloadHandler(
    filename = function() {
      paste(input$text,"_",input$ReadsData, "_Reads.csv", sep = "")
    },
    content = function(file) {
      
      peak.dist.data<-rpeak.dist.data()
 
      if(input$ReadsData=="Glu"){ out<-peak.dist.data$tflist_Glu}
      if(input$ReadsData=="Nit"){ out<-peak.dist.data$tflist_Nit}
      if(input$ReadsData=="Eth"){ out<-peak.dist.data$tflist_Eth}
      if(input$ReadsData=="Ana"){ out<-peak.dist.data$tflist_Ana}
        write.csv(out, file, row.names = FALSE)
    }
  )
  output$downloadPeakData <- downloadHandler(
    filename = function() {
      paste(input$TFs, "_SumPeaks.csv", sep = "")
    },
    content = function(file) {

      out<-read.csv(paste("TF_data_files/GEMPeaks/",input$TFs,"_geneTargetList_190314.csv",sep=""))
      write.csv(out, file, row.names = FALSE)
    }
  )
  output$downloadPeakDataAd <- downloadHandler(
    filename = function() {
      paste(input$TFs, "_All_PeakPos.csv", sep = "")
    },
    content = function(file) {
      out<-read.csv(paste("TF_data_files/GEMPeaks/",input$TFs,"_GEManalysis_190314.csv",sep=""))
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  output$plot1 <- renderPlot({
    peak.dist.data<-rpeak.dist.data()
    TPM<-rTPM()
    BS<-rBS()
    TATA_GTF.df<-rTATA_GTF()
    if (length(peak.dist.data$tflist_Glu)==0)return(NULL)
    p1<-profile_plot_func(rTATA_GTF(),peak.dist.data$tflist_Glu, x,rCurr.Geneseq(),ranges,colorset.TF, BS$d.Glu, rATGC(), input, ryRanges(),rGene_start(),TPM$Glu,rMotiffinder())
    p1
  },bg="transparent")
  output$plot2 <- renderPlot({
    peak.dist.data<-rpeak.dist.data()
    TPM<-rTPM()
    BS<-rBS()
    TATA_GTF.df<-rTATA_GTF()
    if (length(peak.dist.data$tflist_Nit)==0)return(NULL)
    p2<-profile_plot_func(rTATA_GTF(),peak.dist.data$tflist_Nit, x,rCurr.Geneseq(),ranges,colorset.TF, BS$d.Nit, rATGC(), input, ryRanges(),rGene_start(),TPM$Nit,rMotiffinder())
    p2
  },bg="transparent")
  output$plot3 <- renderPlot({
    peak.dist.data<-rpeak.dist.data()
    TPM<-rTPM()
    BS<-rBS()
    if (length(peak.dist.data$tflist_Eth)==0)return(NULL)
    p3<-profile_plot_func(rTATA_GTF(),peak.dist.data$tflist_Eth, x,rCurr.Geneseq(),ranges,colorset.TF, BS$d.Eth, rATGC(), input, ryRanges(),rGene_start(),TPM$Eth,rMotiffinder())
    p3
  },bg="transparent")
  output$plot4 <- renderPlot({
    peak.dist.data<-rpeak.dist.data()
    TPM<-rTPM()
    BS<-rBS()
    if (length(peak.dist.data$tflist_Ana)==0)return(NULL)
    p4<-profile_plot_func(rTATA_GTF(),peak.dist.data$tflist_Ana, x,rCurr.Geneseq(),ranges,colorset.TF, BS$d.Ana, rATGC(), input, ryRanges(),rGene_start(),TPM$Ana,rMotiffinder())
    p4
  },bg="transparent")
  output$plot5 <-renderPlot({
    n=length(name.TF)
    n2=10
    x1<-matrix(0,1,n)
    x2<-matrix(0,1,n)
    y1<-matrix(0,1,n)
    y2<-matrix(0,1,n)
    t<-colorset.TF.df$name.TF
    colval<-colorset.TF.df$myColors
    d=data.frame(t(x1), t(x2), t(y1), t(y2), t, colval)
    colnames(d)<-c("x1","x2","y1","y2","t","colval")
    p1<-ggplot() +
      scale_x_continuous(name="x") +
      scale_y_continuous(name="y")+
      labs(fill = "TF BS", alpha=0.5)+
      scale_fill_manual(values=colorset.TF[1:length(name.TF)])+ 
      geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), alpha=alphaval) +
      theme_void()+
      guides(fill = guide_legend(nrow = nrow(colorset.TF.df), ncol=1))+
      theme(legend.text=element_text(size=14),
            legend.title=element_text(size=16),
            panel.grid.major = element_blank() # get rid of major grid
            , panel.grid.minor = element_blank() # get rid of minor grid
      )
    
    p1
  },bg="transparent")

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
    
    return(paste(paste("Genename Systematic: ", curr.geneSys), paste("Genename Common: ", input.gene),paste("Gene start: ",Chr,ChrPos),paste("Gene up: ", prevgeneC,prevgeneS),paste("Gene down: ", nextgeneC, nextgeneS),sep="\n"))
    }
    else{return("Error: Gene doesn't exist")}
  })
  output$Sequnce_out <- renderText({
    if(input$SeqFind==TRUE){
      curr.geneSys<-gene.Systematic[rGeneidx(),1]
      Geneidx.seq<-grep(paste("^",curr.geneSys,"$",sep=""), as.character(Geneseq[,1]))
      Curr.Geneseq<-Geneseq[Geneidx.seq,-(1:2)]   
      Curr.Geneseq<-as.character(Curr.Geneseq)
      Curr.Geneseq<-sapply(seq(from=1, to=nchar(Curr.Geneseq), by=1), function(i) substr(Curr.Geneseq, i, i))
      Curr.Geneseq<-Curr.Geneseq[1000+(input$SeqMin:input$SeqMax)]
      d<-floor(((1000+as.numeric(input$SeqMax))-(1000+as.numeric(input$SeqMin)))/20)
      Out.Geneseq<-data.frame(matrix(ncol=d,nrow=1))
      o=0
      for (i in 1:d){
        if(i>1){
         o=1}
        
      Out.Geneseq[,i]<-paste(Curr.Geneseq[((i-1)*20+o):(i*20)],collapse="")
      }
      Out.Geneseq<-as.character(Out.Geneseq)
      return(Out.Geneseq)
    }
  })
  rCCM<-reactive({
    val1=0
    if(input$val1=="Exclude dubious and non-Yeast 8" ){val1=1}
    if(input$val1=="Include dubious and non-Yeast 8"){val1=2}
    return(val1)
  })
  output$goterms<-renderDataTable({
    pathways.name<-input$goterm
    pathways.name<-tolower(pathways.name)
    pathways.name<-strsplit(pathways.name, "[+]")
    pathways.name<-unlist(pathways.name)
    pathways.go <- filter(go_bio, grepl(paste(pathways.name, collapse="|"),name_1006))
    pathways.go1<-data.frame(unique(pathways.go$name_1006))
    colnames(pathways.go1)<-"GO-terms"
    return(pathways.go1)
  }, options=list(pageLength=5))
  output$selectInput <- renderText({
    if(input$test=="Fisher"){
      return("The Fisher test is based on the occurens of one transcription factor at the same gene as another transcription factor. 
             The probability is plotted in a heatmap where * denotes p<0.001")
    }
    
    if(input$test=="Heatmap"){
      return("Generates a heatmap of the selected genes")
    }
    if(input$test=="Network"){
      return("The Network shows each transcription factor and its targets vs all other transcription factors. 
             Transcription factors that are highly connected appear closer on the chart. The node size is weighted.")
    }
    if(input$test=="Cluster"){
      return("The genes in the GO-term is clustered by k-medoids clustering based on the transcription factor occurens. Positive medoid coefficient indicates presens of TF while negative coefficient indicates absense of TF in the cluster.")
    }
    if(input$test=="Linear Model"){
      return("Based on the transcription factor occurens a linear model is used to predict the resulting TPM values.")
    }
    })
  plotdata1_p3<-eventReactive(input$search, {
   
    val1<-rCCM()
    val2<-1
    
    pathways.name<-input$goterm
    pathways.name<-tolower(pathways.name)
    pathways.name<-strsplit(pathways.name, "[+]")
    pathways.name<-unlist(pathways.name)
    pathwayData<-goterm_func(pathways.name)
    totalgenes<-nrow(pathwayData)
    ####The data
    conditionname<-"Glu"
    data.list<-data_func(conditionname,pathwayData,val1, val2)
    
    data1<-data.list$data1
    datacondname<-paste(conditionname,".lim",sep="")
    data1<-merge(data1, TPM_data, by="Gene")
    y.ori<-data1[datacondname]
    x.ori<-data1[,name.TF]
    data.list2=data_treatment_func(y.ori,x.ori, data1, val2, conditionname)
    txtstr<-paste( "Selected genes", nrow(data.list2$x),"of total",totalgenes)
    data<-cbind(data.list2$y,data.list2$x)
    ###Fisher
    if(input$test=="Fisher"){
      p1<-tryCatch({fishers_test_func(data.list2$x, txtstr,1)},
                   error=function(cond){
                     return(ggplot()+theme_void()+geom_text(aes(x=1,y=1,label="Not enough genes selected to do statistical function",hjust=(0.4) )))}
      )
    }
    if(input$test=="Heatmap"){
      ##### Only data heatmap
      data<-cbind(data.list2$y,data.list2$x)
      p1<-heatmap_func(data,txtstr,1)

    }
    if(input$test=="Network"){
      p1<-net_func(data.list2$x, txtstr)
      
    }
    if(input$test=="Cluster"){
      p1<-cluster_func(data.list2$x,1, input$slider1)
      
    }
    if(input$test=="Linear Model"){
      datalin<-cbind(data.list2$y,data.list2$x)
      #zero interaction
      p1<- model_zero_func(datalin, txtstr,1)
    }
    
     p1
  })
  plotdata2_p3<-eventReactive(input$search, {
    val1<-rCCM()
    val2<-1
    
    pathways.name<-input$goterm
    pathways.name<-tolower(pathways.name)
    pathways.name<-strsplit(pathways.name, "[+]")
    pathways.name<-unlist(pathways.name)
    pathwayData<-goterm_func(pathways.name)
    totalgenes<-nrow(pathwayData)
    ####The data
    conditionname<-"Nit"
    data.list<-data_func(conditionname,pathwayData,val1, val2)
    data1<-data.list$data1
    datacondname<-paste(conditionname,".lim",sep="")
    data1<-merge(data1, TPM_data, by="Gene")
    y.ori<-data1[datacondname]
    x.ori<-data1[,name.TF]
    
    data.list2=data_treatment_func(y.ori,x.ori, data1, val2, conditionname)
    txtstr<-paste( "Selected genes", nrow(data.list2$x),"of total",totalgenes)
    data<-cbind(data.list2$y,data.list2$x)
    ###Fisher
    if(input$test=="Fisher"){
      p1<-tryCatch({fishers_test_func(data.list2$x, txtstr,1)},
                   error=function(cond){
                     return(ggplot()+theme_void()+geom_text(aes(x=1,y=1,label="Not enough genes selected to do statistical function",hjust=(0.4) )))}
      )
    }
    if(input$test=="Heatmap"){
      data<-cbind(data.list2$y,data.list2$x)
      p1<-heatmap_func(data,txtstr,1)
    }
    if(input$test=="Network"){
      p1<-net_func(data.list2$x, txtstr)
      
    }
    if(input$test=="Cluster"){
      p1<-cluster_func(data.list2$x,1, input$slider1)
    }
    if(input$test=="Linear Model"){
      datalin<-cbind(data.list2$y,data.list2$x)
      #zero interaction
      p1<- model_zero_func(datalin, txtstr,1)
    }
    
    p1
    
  })
  plotdata3_p3<-eventReactive(input$search, {
    val1<-rCCM()
    val2<-1
    
    pathways.name<-input$goterm
    pathways.name<-tolower(pathways.name)
    pathways.name<-strsplit(pathways.name, "[+]")
    pathways.name<-unlist(pathways.name)
    pathwayData<-goterm_func(pathways.name)
    totalgenes<-nrow(pathwayData)
    ####The data
    conditionname<-"Eth"
    data.list<-data_func(conditionname,pathwayData,val1, val2)
    data1<-data.list$data1
    datacondname<-paste(conditionname,".lim",sep="")
    data1<-merge(data1, TPM_data, by="Gene")
    y.ori<-data1[datacondname]
    x.ori<-data1[,name.TF]
    
    data.list2=data_treatment_func(y.ori,x.ori, data1, val2, conditionname)
    txtstr<-paste( "Selected genes", nrow(data.list2$x),"of total",totalgenes)
    data<-cbind(data.list2$y,data.list2$x)
    ###Fisher
    if(input$test=="Fisher"){
      p1<-tryCatch({fishers_test_func(data.list2$x, txtstr,1)},
                   error=function(cond){
                     return(ggplot()+theme_void()+geom_text(aes(x=1,y=1,label="Not enough genes selected to do statistical function",hjust=(0.4) )))}
      )
    }
    if(input$test=="Heatmap"){
      data<-cbind(data.list2$y,data.list2$x)
      p1<-heatmap_func(data,txtstr,1)
      
    }
    if(input$test=="Network"){
      p1<-net_func(data.list2$x, txtstr)
      
    }
    if(input$test=="Cluster"){
      p1<-cluster_func(data.list2$x,1, input$slider1)
    }
    if(input$test=="Linear Model"){
      datalin<-cbind(data.list2$y,data.list2$x)
      #zero interaction
      p1<- model_zero_func(datalin, txtstr,1)
    }
    p1
    
  })
  plotdata4_p3<-eventReactive(input$search, {
    val1<-rCCM()
    val2<-1
    
    pathways.name<-input$goterm
    pathways.name<-tolower(pathways.name)
    pathways.name<-strsplit(pathways.name, "[+]")
    pathways.name<-unlist(pathways.name)
    pathwayData<-goterm_func(pathways.name)
    totalgenes<-nrow(pathwayData)
    ####The data
    conditionname<-"Ana"
    data.list<-data_func(conditionname,pathwayData,val1, val2)
    data1<-data.list$data1
    datacondname<-paste(conditionname,".lim",sep="")
    data1<-merge(data1, TPM_data, by="Gene")
    y.ori<-data1[datacondname]
    x.ori<-data1[,name.TF]
    
    data.list2=data_treatment_func(y.ori,x.ori, data1, val2, conditionname)
    txtstr<-paste( "Selected genes", nrow(data.list2$x),"of total",totalgenes)
    data<-cbind(data.list2$y,data.list2$x)
    ###Fisher
    if(input$test=="Fisher"){
      p1<-tryCatch({fishers_test_func(data.list2$x, txtstr,1)},
                   error=function(cond){
                     return(ggplot()+theme_void()+geom_text(aes(x=1,y=1,label="Not enough genes selected to do statistical function",hjust=(0.4) )))}
      )
    }
    if(input$test=="Heatmap"){
      data<-cbind(data.list2$y,data.list2$x)
      p1<-heatmap_func(data,txtstr,1)
    }
    if(input$test=="Network"){
      p1<-net_func(data.list2$x, txtstr)
      
    }
    if(input$test=="Cluster"){
      p1<-cluster_func(data.list2$x,1, input$slider1)
      
    }
    if(input$test=="Linear Model"){
      datalin<-cbind(data.list2$y,data.list2$x)
      #zero interaction

      p1<- model_zero_func(datalin, txtstr,1)
    }
    p1
    
  })
  
  rdowloadStatData<-reactive({
    val1<-rCCM()
    val2<-1
    
    pathways.name<-input$goterm
    pathways.name<-tolower(pathways.name)
    pathways.name<-strsplit(pathways.name, "[+]")
    pathways.name<-unlist(pathways.name)
    pathwayData<-goterm_func(pathways.name)
    totalgenes<-nrow(pathwayData)
    ####The data
    conditionname<-input$StatData
    data.list<-data_func(conditionname,pathwayData,val1, val2)
    
    data1<-data.list$data1
    datacondname<-paste(conditionname,".lim",sep="")
    data1<-merge(data1, TPM_data, by="Gene")
    y.ori<-data1[datacondname]
    x.ori<-data1[,name.TF]
    data.list2=data_treatment_func(y.ori,x.ori, data1, val2, conditionname)
    txtstr<-paste( "Selected genes", nrow(data.list2$x),"of total",totalgenes)
    data<-cbind(data.list2$y,data.list2$x)

    ###Fisher
    if(input$test=="Fisher"){
      p1<-fishers_test_func(data.list2$x, txtstr,2)
      
      }
    if(input$test=="Heatmap"){
      ##### Only data heatmap
      data<-cbind(data.list2$y,data.list2$x)
      p1<-heatmap_func(data,txtstr,2)
    }
    if(input$test=="Network"){
      p1<-net_func(data.list2$x, txtstr)
    }
    if(input$test=="Cluster"){
      p1<-cluster_func(data.list2$x,2,input$slider1)
    }
    if(input$test=="Linear Model"){
      datalin<-cbind(data.list2$y,data.list2$x)
      #zero interaction
      p1<- model_zero_func(datalin, txtstr,2)
    }
    out.data<-p1
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
  output$plot1_p3 <- renderPlot({
    par(mar=c(1,1,1,1))
    p1<-plotdata1_p3()
    p1
  },bg="transparent")
  output$plot2_p3 <- renderPlot({
    par(mar=c(1,1,1,1))
    p1<-plotdata2_p3()
    p1
  },bg="transparent")
  output$plot3_p3 <- renderPlot({
    par(mar=c(1,1,1,1))
    p3<-  plotdata3_p3()
    p3
  },bg="transparent")
  output$plot4_p3 <- renderPlot({
    par(mar=c(1,1,1,1))
    p4<-    plotdata4_p3()
    p4
  },bg="transparent")
  
  
    }

