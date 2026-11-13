source('09_harvesting_model.R')
# Independent comparison of support-margin signs with the envelope contact interval.
code<-paste(readLines('09_harvesting_core.cpp'),collapse='\n')
diagnostic<-'
// [[Rcpp::export]]
NumericMatrix check_support(NumericMatrix z){
 NumericMatrix v(z.nrow(),4);double B=50,H=25;
 for(int i=0;i<z.nrow();i++){
  double a=z(i,0),b=z(i,1),p=z(i,2);HEnvelope en(a,b,p,B);
  double wB=hw(B,a,b,p),wH=hw(H,a,b,p),d0=hd(0,a,b,p),dH=hd(H,a,b,p),dB=hd(B,a,b,p);
  double A=(b-1)*(1-a),C=a*(b*p-1),c0=A+C,cB=A+C*std::pow(1+B,b);
  double r1=c0>0&&cB<0?wH-H*dH:B*d0-wB;
  double r2=c0<0&&cB>0?(B-H)*dH-(wB-wH):wB-B*dB;
  v(i,0)=r1;v(i,1)=r2;
  v(i,2)=!en.chord&&std::min(en.hi,H)>en.lo+1e-8;
  v(i,3)=!en.chord&&en.hi>std::max(en.lo,H)+1e-8;
 }return v;
}'
env<-new.env();Rcpp::sourceCpp(code=paste(code,diagnostic,sep='\n'),env=env)
grid<-as.matrix(expand.grid(alpha=exp(seq(-3,3,length.out=31)),beta=c(.1,.4,.5,.8,.9),phi=seq(.1,3.1,length.out=61)))
v<-env$check_support(grid)
for(k in 1:2){ok<-abs(v[,k])>1e-7;bad<-which(ok & ((v[,k]>0)!=as.logical(v[,k+2])));if(length(bad)){print(head(cbind(grid[bad,,drop=FALSE],v[bad,,drop=FALSE]),12));stop('support mismatch class ',k,' count ',length(bad))}}
cat('Signed support margins agree with envelope contact intervals at',nrow(grid),'settings.\n')
maxE<-maxr<-0
for(b in c(.1,.5,.9)){
 d<-readRDS(paste0('09_harvesting_results/phi_beta_',b,'.rds'))
 for(phi in list(c(.1,.1),c(.1,3.1),c(3.1,.1),c(3.1,3.1),c(1.6,1.6))){
  i<-which.min(abs(d$axis-phi[1]));j<-which.min(abs(d$axis-phi[2]))
  p<-harvest_population(b,d$axis[i],d$axis[j],order=96L,inner=32L,start=d$ratio_responsive[i,j])
  err<-max(abs(c(p$first[1]-d$E_first[i,j],p$responsive[1]-d$E_responsive[i,j])))/50
  maxE<-max(maxE,err);maxr<-max(maxr,abs(log(p$responsive_ratio/d$ratio_responsive[i,j])))
  cat('phi check',b,phi,'share error',err,'\n')
  stopifnot(err<.001)
 }
 stopifnot(max(d$residual)<1e-7)
}
cat('Largest checked phi-grid time-share error:',maxE,'; log-ratio change:',maxr,'\n')
