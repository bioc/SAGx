# Per Broberg 21FEB01, updated 30JUN05 #
R2mim<-function(probe="12345_at",N = 10,data=inm,out="u:/study/copd/mimscr.txt"){
   lprob<-nchar(probe) # Extract data regarding the candidate probe # 
   x<-data[(substring(paste(data[,1]),1,lprob)==paste(probe) ),];x<-t(x[-1])
   y<-data;y<-t(y[,-1])
   corrs<-cor(x,y) # Calculate the correlation with all probes #
   titt<-cbind(1:nrow(data),data.frame(data[,1]),t(corrs))
   tittut<-titt[sort.list(-abs(titt[,3])),] # Sort the list wrt absolute correlation #
   val<-matrix(tittut[1:N,1]) # Exctract the row numbers of chosen probes #
 # Den här är skum #
   valda<-data[val,]
   x<-t(valda[,-1])
   ut<-cov(x,y=x);utcorrs<-cor(x,y=x)
   mes<-apply(valda[,-1],1,mean)
   v<-ut[upper.tri(ut,diag=TRUE)]
   nz<-dim(y)[1]
   utfil<-c(nz,mes,v)
   ttt<-paste(data.frame(letters)[1:N,1])
   temp<-character(1)
   for (i in 1:N)temp<-paste(temp,ttt[i],sep="")
   com1<-paste("cont ",temp)
   com2<-paste("statread",temp)
   write.table(com1,file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
   probes<-data[val,1]
   ready<-FALSE;start<-1;end<-min(4,N)
   while(!ready){
      if (end==N)ready<-TRUE
          write.table("labels ",file=paste(out),append=TRUE,quote=FALSE,eol=" ",col.names=FALSE,row.names=FALSE)
          for (i in start:end){
             write.table(paste(" ",ttt[i]," "),file=paste(out),quote=FALSE,append=TRUE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ")
             write.table(probes[i],file=paste(out),quote=TRUE,append=TRUE,eol=" ",col.names=FALSE,row.names=FALSE,sep=" ") 
          }
      write.table("",file=paste(out),quote=FALSE,append=TRUE,eol="\n",col.names=FALSE,row.names=FALSE)
      start<-min(start+4,N);end<-min(end+4,N)
   }
   write.table(com2,file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   write.table(utfil,file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   write.table("!",file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   write.table("SatModel;Fit;StepWise",file=paste(out),sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE,append=TRUE)
   colnames(utcorrs)<-paste(probes)
   rownames(utcorrs)<-paste(probes)
   return(utcorrs)
}
