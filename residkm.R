residkm<-function(data,groupcolumn="Cohort",krange=2:10,ksel=T,altfeatnames=NULL,feattypes=NULL,impncp=5,
                impgrouping=NULL,imptypes=NULL,nbcindex="tracew",method=c("kmeans","pam","spectral")){
  library(ggplot2)
  getmode<-function(v){
    unique(v)[which.max(tabulate(match(v,unique(v))))]
  }
  coredat<-data[,-which(names(data)==groupcolumn)]
  colclasses<-sapply(coredat,class)
  if(any(is.na(coredat))){
    naflag<-T
    nadata<-data
    if(is.null(impgrouping)) impgrouping<-rep(1,length=ncol(coredat))
    if(is.null(imptypes)){
      imptypes<-rep("s",ncol(coredat))
      imptypes[which(colclasses%in%c("character","logical","factor"))]<-"n"
    }
    groupnames<-paste0("Group",1:length(impgrouping))
    coredat<-missMDA::imputeMFA(coredat,group=impgrouping,type=imptypes,
                                ncp=impncp,name.group=groupnames,graph=F)$completeObs
    # https://francoishusson.wordpress.com/2017/08/05/can-we-believe-in-the-imputations/
    data<-cbind(coredat,data[,groupcolumn])
    names(data)[ncol(data)]<-groupcolumn
  } else {naflag<-F}
  cohresid<-as.data.frame(matrix(0,nrow(coredat),ncol(coredat)))
  names(cohresid)<-names(coredat)
  for(j in 1:ncol(coredat)){
    form<-formula(paste0(names(coredat)[j],"~",groupcolumn))
    if(!class(data[,names(coredat)[j]])%in%c("numeric","integer")){
      templm<-glm(form,data,family="binomial",na.action=na.exclude)
      tempresid<-resid(templm,type="response")
    } else {
      templm<-lm(form,data,na.action=na.exclude)
      tempresid<-resid(templm)
    }
    cohresid[,names(coredat)[j]]<-tempresid
  };rm(j,templm,tempresid,form)
  if(ksel){
    cohnbc<-NbClust(scale(cohresid),distance = "euclidean",
                    min.nc = min(krange),max.nc = max(krange),
                    method = "complete",index = nbcindex)#$Best.nc[1]
    if(nbcindex=="all"|length(nbcindex)>1){
      bestnc<-getmode(cohnbc$Best.nc[1,])
    } else {
      bestnc<-cohnbc$Best.nc[1]
    }
  } else {
    cohnbc<-krange[1]
    bestnc<-krange[1]
  }
  if(length(method)>1) method<-method[1]
  if(method=="kmeans"){
    km<-kmeans(scale(cohresid),bestnc,iter.max=100,nstart=1000)
  } else if(method=="pam"){
    km<-fpc::pamk(scale(cohresid),krange=bestnc,criterion="ch")
    km$cluster<-km$pamobject$clustering
    km$centers<-km$pamobject$medoids
    rownames(km$centers)<-1:nrow(km$centers)
  } else if(method=="specc"){
    km<-fpc::speccCBI(scale(cohresid),bestnc)
    km$cluster<-km$partition
  }
  
  reorder_clusters<-function(x){
    fx<-as.factor(x)
    summfx<-summary(fx)
    newlabs<-summfx[order(-summfx)]
    #levels(fx)<-newlabs#names(newlabs)
    fx<-factor(fx,levels=names(newlabs))
    return(as.numeric(fx))
  }
  origclustcounts<-summary(as.factor(km$cluster))
  km$cluster<-reorder_clusters(km$cluster)
  newclustcounts<-summary(as.factor(km$cluster))
  km$centers<-km$centers[match(newclustcounts,origclustcounts),]
  rownames(km$centers)<-1:nrow(km$centers)
  kmc<-as.data.frame(km$centers)
  if(!is.null(altfeatnames)){
    featnames<-altfeatnames
  } else {
    featnames<-names(coredat)
  }
  names(kmc)<-featnames
  kmc$Cluster<-as.factor(paste0("Cluster",rownames(km$centers)))
  clusternobs<-summary(as.factor(km$cluster))
  levels(kmc$Cluster)<-paste0(levels(kmc$Cluster)," (",clusternobs,")")
  kmc<-reshape2::melt(kmc,id.vars="Cluster")
  kmc$variable<-factor(kmc$variable,levels=rev(featnames))
  if(!is.null(feattypes)){
    ftd<-data.frame(Types=feattypes,Feats=featnames)
    kmc$Type<-factor(ftd[match(as.character(kmc$variable),ftd$Feats),"Types"],levels=unique(feattypes))
    centerplot<-ggplot(kmc,aes(y=variable,x=value))+theme_bw()+
      geom_bar(aes(fill=value),stat="identity",alpha=0.7)+
      geom_vline(xintercept=0)+geom_text(aes(label=round(value,2)))+
      scale_fill_viridis_c()+facet_nested(Cluster+Type~.,scales="free_y",space="free_y")+
      xlab("Center Z-score")+
      theme(legend.position="none",axis.title.y=element_blank())
  } else {
    centerplot<-ggplot(kmc,aes(y=variable,x=value))+theme_bw()+
      geom_bar(aes(fill=value),stat="identity",alpha=0.7)+
      geom_vline(xintercept=0)+geom_text(aes(label=round(value,2)))+
      scale_fill_viridis_c()+facet_grid(Cluster~.)+xlab("Center Z-score")+
      theme(legend.position="none",axis.title.y=element_blank())
  }
  
  summarybycluster<-function(data,clustvec,contvars,catvars,
                             contwhitecutoff=0,catwhitecutoff=0,scale=T,
                             title=NULL,titlesize=12,xlabelsize=12,
                             yaxistextsize=10,xaxistextsize=10,na.rm=T){
    data$Cluster<-as.factor(clustvec)
    countsumm<-summary(data$Cluster)
    clustcountdf<-data.frame(Cluster=names(countsumm),Count=countsumm)
    data$Count<-clustcountdf[match(data$Cluster,clustcountdf$Cluster),"Count"]
    data$Cluster<-as.factor(paste0(as.character(data$Cluster),
                                   " (n=",data$Count,")"))
    data$Count<-NULL
    datacat_long<-reshape2::melt(data[,which(!names(data)%in%contvars)],
                                 id.vars="Cluster")
    heatdfcat<-aggregate(
      value~Cluster+variable,datacat_long,
      function(x) (length(which(x%in%c("Y",1)))/length(x))*100,na.action=NULL)
    datacont<-data[,c(contvars)]
    if(any(scale)){
      if(length(scale)==1){scale<-rep(scale,length(contvars))}
      for(i in 1:length(contvars)){
        if(scale[i]){datacont[,contvars[i]]<-as.numeric(scale(datacont[,contvars[i]]))} 
      }
    }
    datacont$Cluster<-data$Cluster
    #datacont[,1:5]<-sweep(datacont[,1:5],2,contcenters,FUN="+")
    datacont_long<-reshape2::melt(datacont,id.vars="Cluster")
    heatdfcont<-aggregate(value~Cluster+variable,datacont_long,mean)
    heatdfcont2<-aggregate(value~Cluster+variable,datacont_long,sd)
    heatdfcont$SD<-heatdfcont2$value
    rm(heatdfcont2)
    heatdfcont$TextLabel<-paste0(round(heatdfcont$value,1),
                                 " (",round(heatdfcont$SD,1),")")
    heatdfcat$variable<-factor(heatdfcat$variable,
                               levels=rev(levels(heatdfcat$variable)))
    heatdfcont$variable<-factor(heatdfcont$variable,
                                levels=rev(levels(heatdfcont$variable)))
    heatdfcont$TextColor<-ifelse(heatdfcont$value>contwhitecutoff,"black","white")
    heatdfcat$TextColor<-ifelse(heatdfcat$value>catwhitecutoff,"black","white")
    if(any(scale)){
      heatdfcont$StripTitle<-"Mean (SD) Z Score"
    } else {
      heatdfcont$StripTitle<-"Mean (SD)"
    }
    heatdfcat$StripTitle<-"Percent"
    
    if(na.rm){
      heatdfcont<-heatdfcont[heatdfcont$Cluster!="NA",]
      heatdfcat<-heatdfcat[heatdfcat$Cluster!="NA",]
    }
    
    g1<-ggplot(heatdfcont,aes(y=variable,x=Cluster,fill=value))+
      geom_tile(alpha=0.6)+
      geom_text(aes(label=TextLabel,color=TextColor),size=3.5)+
      facet_grid(StripTitle~.,space="free",scales="free")+
      scale_fill_viridis_c()+cowplot::theme_cowplot()+
      scale_color_manual(values=c("black","white"),guide=F)+
      theme(legend.position="none",
            axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text=element_text(size=12),
            axis.text.y=element_text(size=yaxistextsize))
    if(!is.null(title)){g1<-g1+ggtitle(title)+
      theme(plot.title=element_text(size=titlesize,hjust=0.5))}
    g2<-ggplot(heatdfcat,aes(y=variable,x=Cluster,fill=value))+
      geom_tile(alpha=0.6)+
      geom_text(aes(label=round(value,1),color=TextColor),size=4)+
      facet_grid(StripTitle~.,space="free",scales="free")+
      scale_fill_viridis_c(option="C")+cowplot::theme_cowplot()+
      scale_color_manual(values=c("black","white"),guide=F)+
      theme(legend.position="none",
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=xlabelsize),
            strip.text=element_text(size=12),
            axis.text.x=element_text(size=xaxistextsize),
            axis.text.y=element_text(size=yaxistextsize))
    catlength<-length(unique(heatdfcat$variable))
    contlength<-length(unique(heatdfcont$variable))
    catplotht<-(catlength)/(catlength+contlength)
    contplotht<-contlength/(catlength+contlength)
    fig<-cowplot::plot_grid(g1,g2,nrow=2,align="v",
                            rel_heights=c(contplotht,catplotht))
    return(fig)
  }
  contvars<-featnames[which(colclasses%in%c("numeric","integer"))]
  if(exists("nadata")) {
    meanfigdata<-nadata[,-which(names(nadata)==groupcolumn)]
  } else {
    meanfigdata<-data[,-which(names(data)==groupcolumn)]
  }
  names(meanfigdata)<-featnames
  meanfig<-summarybycluster(meanfigdata,km$cluster,
                            contvars,setdiff(featnames,contvars))
  eudist<-function(x) {
    out<-matrix(0,ncol(x),ncol(x))
    for(i in 1:ncol(x)){
      for(j in 1:ncol(x)){
        out[i,j]<-sqrt(sum((x[,i]-x[,j])^2))
      };rm(j)
    };rm(i)
    rownames(out)<-colnames(out)<-colnames(x)
    return(out)
  }
  
  allcenters<-km$centers
  rownames(allcenters)<-paste0("C",1:nrow(km$centers))
  allcenters<-t(allcenters)
  eudists<-eudist(allcenters)
  eudlabs<-apply(eudists,2,function(x) round(x,2))
  hm1<-gplots::heatmap.2(eudists,scale="none",col=hcl.colors(50),trace="none",notecex=1.25,
                         cellnote=eudlabs,notecol="black",tracecol="black")
  
  if(naflag){
    return(list(Kmeans=km,KChoice=cohnbc,ResidualData=cohresid,CenterPlot=centerplot,MeanPlot=meanfig,CenterEuDist=hm1,ImpData=data))
  } else {
    return(list(Kmeans=km,KChoice=cohnbc,ResidualData=cohresid,CenterPlot=centerplot,MeanPlot=meanfig,CenterEuDist=hm1))
  }
}
