source('10_campaign_precision_maps.R')
verify_campaign_precision_maps<-function(){
 for(e in c(.05,.1,.2)){
  d<-assemble_campaign_precision_fields(e);stopifnot(length(d$axis)==25L)
  for(mode in names(d$panels))for(c in 1:3)for(r in 1:3){
   p<-d$panels[[mode]][[c,r]];q<-d$panels[[mode]][[r,c]]
   stopifnot(all(is.finite(p$C1)),all(is.finite(p$C2)),min(p$C1,p$C2)>=-1e-10,
    max(p$C1+p$C2)<=1+1e-9,max(abs(p$C1-t(q$C2)))<1e-12)
   if(c==r)stopifnot(max(abs(diag(p$C1)-diag(p$C2)))<1e-12)
   fine<-campaign_bilinear(p$axis,p$C1+p$C2,seq(.1,.9,length.out=73))
   stopifnot(min(fine)>=min(p$C1+p$C2)-1e-10,max(fine)<=max(p$C1+p$C2)+1e-10)
  }
 }
 message('All parameter maps have valid budgets, exact exchange symmetry and bounded interpolation.')
 invisible(TRUE)
}
if(sys.nframe()==0L)verify_campaign_precision_maps()
