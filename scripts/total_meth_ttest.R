# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE) {
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

# whole-genome
{
  mean_wt = list(CG = 25.8222184233241,
                 CHG = 6.837569844,
                 CHH = 2.264188479)
  
  mean_mto1 = list(CG = 25.94971156,
                   CHG = 8.239533867,
                   CHH = 2.540204156)
  
  
  
  sd_wt = list(CG = 0.067418241,
               CHG = 0.335858464,
               CHH = 0.160429206)
  
  sd_mto1 = list(CG = 0.203546496,
                 CHG = 0.61711384,
                 CHH = 0.119720471)
}

# Heterochromatin
{
  mean_wt = list(CG = 60.47926345,
                 CHG = 21.67824442,
                 CHH = 5.29480464628366)
  
  mean_mto1 = list(CG = 60.7700085259239,
                   CHG = 26.7109355410493,
                   CHH = 6.22095974639316)
  
  
  
  sd_wt = list(CG = 0.170526611,
               CHG = 1.013817966,
               CHH = 0.326478877)
  
  sd_mto1 = list(CG = 0.509630048,
                 CHG = 2.207031562,
                 CHH = 0.363140872)
}

# Euchromatin
{
  mean_wt = list(CG = 14.11139276,
                 CHG = 2.241731896,
                 CHH = 1.222088897)
  
  mean_mto1 = list(CG = 14.17859944,
                   CHG = 2.516469179,
                   CHH = 1.274194071)

  
  
  sd_wt = list(CG = 0.036019641,
               CHG = 0.126253433,
               CHH = 0.103051467)
  
  sd_mto1 = list(CG = 0.105349122,
                 CHG = 0.12407056,
                 CHH = 0.035314212)
}

for (context in c("CG","CHG","CHH")) {
  tt = t.test2(mean_wt[[context]], mean_mto1[[context]],
               sd_wt[[context]], sd_mto1[[context]],
               2, 3, equal.variance = F)
  
  message("pValue in ", context,": ", round(as.numeric(tt["p-value"]),3))

}
################################################
##### old values

# # whole-genome
# {
#   mean_wt = list(CG = 25.8222184233241,
#                  CHG = 6.837569844,
#                  CHH = 2.264188479)
#   
#   mean_mto1 = list(CG = 25.94971156,
#                    CHG = 8.239533867,
#                    CHH = 2.540204156)
#   
#   
#   
#   sd_wt = list(CG = 0.067418241,
#                CHG = 0.335858464,
#                CHH = 0.160429206)
#   
#   sd_mto1 = list(CG = 0.203546496,
#                  CHG = 0.61711384,
#                  CHH = 0.119720471)
# }
# 
# # Heterochromatin
# {
#   mean_wt = list(CG = 64.59572273,
#                  CHG = 23.36609872,
#                  CHH = 5.631386549)
#   
#   mean_mto1 = list(CG = 64.88534458,
#                    CHG = 28.82821874,
#                    CHH = 6.634528996)
#   
#   
#   
#   sd_wt = list(CG = 0.176684015,
#                CHG = 1.095036704,
#                CHH = 0.345368844)
#   
#   sd_mto1 = list(CG = 0.544601226,
#                  CHG = 2.404550818,
#                  CHH = 0.39477986)
# }
# 
# # Euchromatin
# {
#   mean_wt = list(CG = 14.7148213,
#                  CHG = 2.458310784,
#                  CHH = 1.282856463)
#   
#   mean_mto1 = list(CG = 14.78987639,
#                    CHG = 2.78077083,
#                    CHH = 1.34633032)
# 
#   
#   
#   sd_wt = list(CG = 0.040790156,
#                CHG = 0.135553252,
#                CHH = 0.106438259)
#   
#   sd_mto1 = list(CG = 0.111292709,
#                  CHG = 0.142819606,
#                  CHH = 0.03895111)
# }
# 
# for (context in c("CG","CHG","CHH")) {
#   tt = t.test2(mean_wt[[context]], mean_mto1[[context]],
#                sd_wt[[context]], sd_mto1[[context]],
#                2, 3, equal.variance = F)
#   
#   message("pValue in ", context,": ", round(as.numeric(tt["p-value"]),3))
# 
# }
# 
# 