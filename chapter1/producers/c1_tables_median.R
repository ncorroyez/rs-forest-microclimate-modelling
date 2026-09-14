suppressPackageStartupMessages(library(data.table))
set.seed(1);B<-2000
agg<-function(x,fn){x<-x[is.finite(x)];m<-fn(x);b<-replicate(B,fn(sample(x,length(x),TRUE)));sprintf("%+.3f (%+.3f, %+.3f)",m,quantile(b,.025),quantile(b,.975))}
tab<-function(csv,label){D<-fread(csv);D[,P:=factor(P,levels=paste0("P",1:4))]
 rows<-list(c("LAI +0.5","LAI_up",0.5),c("LAI -0.5","LAI_dn",0.5),c("fCover+10","fCov_up",0.10),
            c("fCover-10","fCov_dn",0.10),c("Hmax +1","Hmax_up",1),c("Hmax -1","Hmax_dn",1),c("profile","dT_LAD",1))
 for(fn in c("mean","median")){cat(sprintf("\n===== %s : %s =====\n",label,toupper(fn)))
  FN<-get(fn)
  for(r in rows){cat(sprintf("%-10s",r[1]))
   for(p in paste0("P",1:4)){x<-as.numeric(D[[r[2]]][D$P==p])*as.numeric(r[3]);cat("  ",agg(x,FN))}
   cat("\n")}}}
tab("out_files/Chapter1/tables/perturb_chs41_nowind.csv","TABLE 2 (native ΔTmax)")
sl<-"out_files/Chapter1/tables/perturb_chs41_nowind_SLOPE.csv"
if(file.exists(sl)) tab(sl,"TABLE E1 (slope)") else cat("\n[SLOPE csv missing — skip E1]\n")
