# Additional A-D columns for Supplement Figures S14 and S15.
# An optional alim changes that column's displayed talent range only.
asymmetric_supplement_examples <- list(
  # Orange crosses.
  S14=list(
    list(beta1=.92,beta2=.4,phi1=1.5,phi2=1),
    list(beta1=.8,beta2=.5,phi1=1.5,phi2=1),
    list(beta1=.5,beta2=.9,phi1=1.4,phi2=.3),
    list(beta1=.9,beta2=.9,phi1=2.1,phi2=.5)),
  # Green crosses. C and D show the same population at different zooms.
  S15=list(
    list(beta1=.9,beta2=.5,phi1=2,phi2=1.5),
    list(beta1=.5,beta2=.9,phi1=1.3,phi2=2.2),
    list(beta1=.5,beta2=.1,phi1=1.8,phi2=2,alim=c(0,.5)),
    list(beta1=.5,beta2=.1,phi1=1.8,phi2=2,alim=c(0,10000)))
)
