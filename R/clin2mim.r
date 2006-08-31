# Copyright Per Broberg #
clin2mim<-function(variable="FEV1.ACTUAL",data=dbs,clindat=clinical,probes=probes,N=10,out="mimscr.txt"){
   data<-data.frame(data)
   if(!is.data.frame(clindat)) stop("The clinical data is not a dataframe.")
   data[,1]<-apply(as.data.frame(data[,1]),1,toupper)
   clindat[,1]<-apply(as.data.frame(clindat[,1]),1,toupper)
   clindat<-clindat[,c(1,which(colnames(clindat)==variable))]
   totdata<-merge(clindat,data,by.x=1, by.y=1)
   x<-totdata[,2]
   y<-totdata[,-c(1,2)] # y now holds the probes' expression in columns #
   corrs<-cor(x,y,use="pairwise.complete.obs") # Calculate the correlation with all probes #
   titt<-cbind(1:ncol(y),probes,t(corrs))
   tittut<-titt[sort.list(-abs(titt[,3])),] # Sort the list wrt absolute correlation #
   val<-matrix(tittut[1:(N-1),1]) # Exctract the row numbers of chosen probes #
   z<-t(y);valda<-z[val,]
   w<-t(rbind(x,valda))
   ut<-cov(w,y=w,use="pairwise.complete.obs");utcorrs<-cor(w,y=w,use="pairwise.complete.obs")
   mes<-apply(w,2,function(x) mean(na.omit(x)))
   v<-ut[upper.tri(ut,diag=TRUE)]
   nz<-dim(y)[1]
   utfil<-c(nz,mes,v)
   ttt<-paste(data.frame(letters)[1:N,1])
   temp<-character(1)
   for (i in 1:N)temp<-paste(temp,ttt[i],sep="")
   com1<-paste("cont ",temp)
   com2<-paste("statread",temp)
   write.table(com1,file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
   labs<-c(variable,probes[val])
   ready<-F;start<-1;end<-min(4,N)
   while(!ready){
      if (end==N)ready<-TRUE
          write.table("labels ",file=paste(out),append=TRUE,quote=FALSE,eol=" ",col.names=FALSE,row.names=FALSE)
          for (i in start:end){
             write.table(paste(" ",ttt[i]," "),file=paste(out),quote=FALSE,append=TRUE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ")
             write.table(labs[i],file=paste(out),quote=TRUE,append=TRUE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ") 
          }
      write.table("",file=paste(out),quote=FALSE,append=TRUE,eol="\n",col.names=FALSE,row.names=FALSE)
      start<-min(start+4,N);end<-min(end+4,N)
   }
   write.table(com2,file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   write.table(utfil,file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   write.table("!",file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   write.table("SatModel;Fit;StepWise",file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   colnames(utcorrs)<-paste(labs)
   rownames(utcorrs)<-paste(labs)
   return(utcorrs)
}
