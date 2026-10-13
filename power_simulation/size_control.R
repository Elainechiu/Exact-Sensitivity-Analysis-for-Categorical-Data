##############################################################################################
#                                                                                            #
#   This file discusses the size control, type I error rate from the exact                   #
#   method developed in Section 4, and from the normal approximation method                  #
#   developed in Section 5 under I by two table sign-score tests.                            #
#                                                                                            #
#   The simulation results are presented                                                     #
#   in the paper "Exact Sensitivity Analysis for Observational Studies of Contingency        #
#   Tables" -- Section 7, "Size Control and Power of a Sensitivity Analysis" and             #
#   some additional results are in Appendix F                                                #
#   We include two equivalent forms of sampling from the                                     #
#   distribution under the null hypothesis of no treatment effect but under                  #
#   the worst-case unmeasured confounder under the generic bias sensitivity model            #
#   The first function is more similar to how the treatment model is written;                #
#   the second function is more computationally efficient, applying our Corollary            #
#   1, which shows that the cell counts from an I by two table under the worst-case          #
#   unmeasured confounder of the sign-score test follows a multivariate extended             #
#   hypergeometric distribution, which can be sampled from the R package BiasedUrn.          #
#   Under each data generating process, we compute the upper bound on the p-value            #
#   with the exact method, or compute the upper bound on the p-value with the normal         #
#   approximation, which is to first compute the first two moments and use our Corollary     #
#   2 in Section 5.                                                                          #                                           
#                                                                                            #
#                                                                                            #  
##############################################################################################






########################  Start of the codes ###################################################
################################################################################################


######################
# load the libraries #
######################




library(ggplot2)
library(Rcpp)
source("sensitivityIxJ.R")
library(BiasedUrn)
sourceCpp("sensitivityIxJ.cpp")
library(multicool)















#######################################################################################
#                                                                                     #
# Develop two equivalent forms of sampling sign-score test from I by two tables       #
# under null hypothesis of no treatment effect and generic bias treatment assignment  #
# at the worst-case unmeasured confounder. We further check that these two functions  #
# are truly equivalent.                                                               #
#######################################################################################


### sampling function: form 1 

biased.ordinal.test.multiple.treatment.binary.outcome.moment = function(treatment.margins, outcome.margins, gamma, delta, treatment.scores, outcome.scores, mc.iteration = 10000){
  
  # Check that we have at least 2 treatments
  if(length(treatment.margins) < 2){
    stop("this function requires at least two treatments")
  }
  
  # Check that outcome is binary
  if(length(outcome.margins) != 2){
    stop("this function only discusses binary outcome")
  }
  
  # Check that treatment.scores matches the number of treatments
  if(length(treatment.scores) != length(treatment.margins)){
    stop("treatment.scores must have the same length as treatment.margins")
  }
  
  # Check that delta matches the number of treatments
  if(length(delta) != length(treatment.margins)){
    stop("delta must have the same length as treatment.margins")
  }
  
  ### Create treatment assignment vector
  z.seed = c()
  for(i in 1:length(treatment.margins)){
    z.seed = c(z.seed, rep(i, treatment.margins[i]))
  }
  
  # Generate all permutations
  mc.Obj = initMC(z.seed)
  all.permutation = allPerm(mc.Obj)
  each.z.probability.weight = rep(NA, nrow(all.permutation))
  
  # Worst case outcome allocation
  worst.case.u = c(rep(0, outcome.margins[1]), rep(1, outcome.margins[2]))
  
  # Compute weights for each permutation
  for(row in 1:nrow(all.permutation)){
    ## compute based on the model and the u.allocation, which we use the outcome margins 
    permute.treatment = all.permutation[row,]
    weight = exp(gamma * delta[permute.treatment] %*% worst.case.u)
    each.z.probability.weight[row] = weight
  }
  
  ## here finishes the construction of the weight for each z 
  
  # Create outcome seed
  outcome.seed = c(rep(1, outcome.margins[1]), rep(2, outcome.margins[2]))
  
  # Transform function to compute test statistic
  transform.fun = function(tb){
    scores = outer(treatment.scores, outcome.scores)
    return(sum(tb * scores))
  }
  
  # Monte Carlo simulation
  obs.stat.vec = rep(NA, mc.iteration)
  for(mc in 1:mc.iteration){
    ## resample
    sampled.z.index = sample(x = 1:nrow(all.permutation), size = 1, replace = FALSE, prob = each.z.probability.weight)
    ## 
    sampled.z = all.permutation[sampled.z.index,]
    sampled.tb = table(sampled.z, outcome.seed)
    obs.stat.vec[mc] = transform.fun(sampled.tb)
  }
  
  return(list(obs.stat.vec = obs.stat.vec, obs.mean = mean(obs.stat.vec), obs.sd = sd(obs.stat.vec)))
}



## sampling function: form 2


biased.hyper.ordinal.test.multiple.treatment.binary.outcome.moment = function(
    treatment.margins, 
    outcome.margins, 
    delta, 
    gamma, 
    treatment.scores, 
    outcome.scores, 
    mc.iteration = 10000
){
  # Check that we have at least 2 treatments
  if(length(treatment.margins) < 2){
    stop("this function requires at least two treatments")
  }
  
  # Check that outcome is binary
  if(length(outcome.margins) != 2){
    stop("this function only discusses binary outcome")
  }
  
  # Check that delta has the same length as treatment.margins
  if(length(delta) != length(treatment.margins)){
    stop("delta should have the same length as treatment.margins")
  }
  
  # Check that treatment.scores has the same length as treatment.margins
  if(length(treatment.scores) != length(treatment.margins)){
    stop("treatment.scores should have the same length as treatment.margins")
  }
  
  # Check that delta is nondecreasing
  if(any(diff(delta) < 0)){
    stop("delta should be nondecreasing")
  }
  
  # Check that treatment.scores is nondecreasing
  if(any(diff(treatment.scores) < 0)){
    stop("treatment.scores should be nondecreasing")
  }
  
  transform.fun = function(tb){
    scores = outer(treatment.scores, outcome.scores)
    return(sum(tb * scores))
  }
  
  weights = exp(gamma * delta)
  hyper.obs.stat.vec = rep(NA, mc.iteration)
  
  for(mc in 1:mc.iteration){
    # Sample allocation of outcome 2 across treatments:
    trt.outcome.2 <- as.vector(rMFNCHypergeo(
      nran = 1,
      m = treatment.margins,
      n = outcome.margins[2], # total outcome=2
      odds = weights
    ))
    # Outcome 1 by subtraction:
    trt.outcome.1 <- treatment.margins - trt.outcome.2
    obs.table <- cbind(trt.outcome.1, trt.outcome.2)
    ## compute the test statistic
    hyper.obs.stat.vec[mc] = transform.fun(obs.table)
  }
  
  return(list(obs.stat.vec = hyper.obs.stat.vec, obs.mean = mean(hyper.obs.stat.vec), obs.sd = sd(hyper.obs.stat.vec)))
}


## Do experiment to verify the moments of ordinal test from two function will be almost identical given the same biased treatment assignment 


biased.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,2,3), outcome.margins = c(3,4), delta = c(0,0,1),gamma=0.5,treatment.scores = c(0,1,2), outcome.scores = c(0,1), mc.iteration = 100000)
biased.hyper.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,2,3), outcome.margins = c(3,4), delta = c(0,0,1),gamma=0.5,treatment.scores = c(0,1,2), outcome.scores = c(0,1), mc.iteration = 100000)




biased.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,3,3), outcome.margins = c(4,4), delta = c(0,1,1),gamma=1,treatment.scores = c(0,2,3), outcome.scores = c(0,1), mc.iteration = 100000)
biased.hyper.ordinal.test.multiple.treatment.binary.outcome.moment(treatment.margins = c(2,3,3), outcome.margins = c(4,4), delta = c(0,1,1),gamma=1,treatment.scores = c(0,2,3), outcome.scores = c(0,1), mc.iteration = 100000)






##### The following function samples the signed score test from the biased distribution 
##### for each Monte Carlo simulations, and under each Monte Carlo simulation, 
##### obtain the exact p-value and the normally-approximated p value.
##### Finally return a graph of the type I error rate 

normal.exact.multiple.treatment.binary.outcome.size.simu.fun = function(
    treatment.margins, 
    outcome.margins, 
    delta, 
    gamma, 
    treatment.scores, 
    outcome.scores, 
    mc.iteration = 10000,
    verbose = FALSE,
    seed = 100,
    base_size = 20
){
  set.seed(seed)
  
  # Check that we have at least 2 treatments
  if(length(treatment.margins) < 2){
    stop("this function requires at least two treatments")
  }
  
  # Check that outcome is binary
  if(length(outcome.margins) != 2){
    stop("this function only discusses binary outcome")
  }
  
  # Check that delta has the same length as treatment.margins
  if(length(delta) != length(treatment.margins)){
    stop("delta should have the same length as treatment.margins")
  }
  
  # Check that treatment.scores has the same length as treatment.margins
  if(length(treatment.scores) != length(treatment.margins)){
    stop("treatment.scores should have the same length as treatment.margins")
  }
  
  # Check that delta is nondecreasing
  if(any(diff(delta) < 0)){
    stop("delta should be nondecreasing")
  }
  
  # Check that treatment.scores is nondecreasing
  if(any(diff(treatment.scores) < 0)){
    stop("treatment.scores should be nondecreasing")
  }
  
  if(sum(treatment.margins) != sum(outcome.margins)){
    stop("what is the correct sample size?")
  }
  
  transform.fun = function(tb){
    scores = outer(treatment.scores, outcome.scores)
    return(sum(tb * scores))
  }
  
  weights = exp(gamma * delta)
  hyper.obs.stat.vec = rep(NA, mc.iteration)
  
  
  for(mc in 1:mc.iteration){
    
    # Sample allocation of outcome 2 across treatments:
    trt.outcome.2 <- as.vector(rMFNCHypergeo(
      nran = 1,
      m = treatment.margins,
      n = outcome.margins[2], # total outcome=2
      odds = weights
    ))
    
    # Outcome 1 by subtraction:
    trt.outcome.1 <- treatment.margins - trt.outcome.2
    
    obs.table <- cbind(trt.outcome.1, trt.outcome.2)
    
    ## compute the test statistic
    hyper.obs.stat.vec[mc] = transform.fun(obs.table)
    
  }
  ## after running the first set of monte carlo, we obtain an estimate of the worst.case.mean and the worst.case.sd
  worst.case.mean = mean(hyper.obs.stat.vec)
  worst.case.sd   =   sd(hyper.obs.stat.vec)
  
  ## now starts another round of monte carlo simulation, but only for pvalue
  exact.p.vec = rep(NA, mc.iteration)
  normal.p.vec= rep(NA, mc.iteration)
  ## renew the hyper.obs.stat.vec
  hyper.obs.stat.vec=rep(NA,mc.iteration)
  
  
  for(mc in 1:mc.iteration){
    ## use the same sampling steps 
    trt.outcome.2 <- as.vector(rMFNCHypergeo(
      nran = 1,
      m = treatment.margins,
      n = outcome.margins[2], # total outcome=2
      odds = weights
    ))
    
    trt.outcome.1 <- treatment.margins - trt.outcome.2
    
    obs.table <- cbind(trt.outcome.1, trt.outcome.2)
    
    obs.stat  = transform.fun(obs.table)
    ## compute the test statistic
    hyper.obs.stat.vec[mc] = obs.stat 
    
    ## compute the exact p value 
    exact.p.vec[mc] = exact.score.sen.IxJ(obs.table = obs.table, gamma = gamma, delta = delta, treatment.scores = treatment.scores, outcome.scores = outcome.scores)$max.prob
    
    
    ## make sure that both the exact and the normal uses the same worst-case unmeasured confounder 
    normal.p.vec[mc] = norm.score.sen.IxJ(obs.table = obs.table, gamma = gamma, delta = delta, treatment.scores = treatment.scores, outcome.scores = outcome.scores, u_space = t(as.matrix(c(0, outcome.margins[2]))))$max.prob
    
    if(mc %% 10 == 0){
      cat("current mc", mc, "\n")
    }
    
    
  }
  
  ## after all the running, return
  test.stat.hist = hist(hyper.obs.stat.vec)
  
  ## construct the type I error rate graph 
  alpha_seq <- seq(0, 1, length.out = mc.iteration)
  df_exact <- data.frame(
    alpha = alpha_seq,
    type1err = sapply(alpha_seq, function(a) mean(exact.p.vec <= a)),
    Method = "Exact"
  )
  df_normal <- data.frame(
    alpha = alpha_seq,
    type1err = sapply(alpha_seq, function(a) mean(normal.p.vec <= a)),
    Method = "Normal"
  )
  
  df_uniform <- data.frame(
    alpha = alpha_seq,
    type1err = alpha_seq,
    Method = "Uniform"
  )
  
  ## report the max of the inflation of type I error rate, which is type1err - alpha_seq
  normal.type.1.err = sapply(alpha_seq, function(a) mean(normal.p.vec <=a))
  
  # Combine data frames
  plot_df <- rbind(df_exact, df_normal, df_uniform)
  plot_df$Method <- factor(plot_df$Method, levels = c("Exact", "Normal", "Uniform"))
  
  # Plotting
  p <- ggplot(plot_df, aes(x = alpha, y = type1err, color = Method, linetype = Method)) +
    geom_line(size = 1.4) +
    scale_color_manual(values = c("#F8766D", "#00BFC4", "black")) +
    scale_linetype_manual(values = c("solid", "longdash", "dashed")) +
    labs(x = "Nominal Significance Level", y = "Empirical Rejection Rate") +
    theme_minimal(base_size = base_size, base_family = "Times New Roman") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = base_size, family = "Times New Roman"),
      axis.title = element_text(size = base_size, family = "Times New Roman"),
      axis.text = element_text(size = base_size, family = "Times New Roman"),
      legend.key.width = unit(2.5, "cm"),
      legend.position = "bottom"
    )
  return(list(obs.stat.vec = hyper.obs.stat.vec, test.stat.hist = test.stat.hist, type.1.error.graph = p, worst.case.mean = worst.case.mean, worst.case.sd = worst.case.sd, exact.p.vec = exact.p.vec, normal.p.vec = normal.p.vec))
  
  
}





######### A helper function to combine plots 

library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)

combine_plots_patchwork <- function(plot1, plot2, plot3,
                                    titles = c("Scenario 1", "Scenario 2", "Scenario 3"),
                                    base_size = 18) {  # Adjust this to match your axis text
  
  # Modify plots
  p1 <- plot1 + 
    theme(axis.title.x = element_blank()) +
    ggtitle(titles[1]) +
    theme(plot.title = element_text(hjust = 0.5, size = base_size)) +
    coord_fixed(ratio = 1)  # Make plot square
  
  p2 <- plot2 + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle(titles[2]) +
    theme(plot.title = element_text(hjust = 0.5, size = base_size)) +
    coord_fixed(ratio = 1)  # Make plot square
  
  p3 <- plot3 + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle(titles[3]) +
    theme(plot.title = element_text(hjust = 0.5, size = base_size)) +
    coord_fixed(ratio = 1)  # Make plot square
  
  # Combine with patchwork and add shared x-axis label with matching font
  combined <- (p1 | p2 | p3) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = "Nominal Significance Level",
      theme = theme(plot.caption = element_text(hjust = 0.5, size = base_size, family = "Times New Roman"))
    ) &
    theme(legend.position = "right",
          legend.text = element_text(size = base_size),    # Match axis text size
          legend.title = element_text(size = base_size),   # Match axis text size
          legend.key.size = unit(0.3, "cm"),              # Smaller legend keys
          legend.key.width = unit(1.4, "cm"),             # Shorter legend lines
          legend.margin = margin(l = -10))                 # Reduce legend margin
  
  return(combined)
}




combined_plot_clean <- combine_plots_patchwork(
  size_control_res_1$type.1.error.graph,
  size_control_res_2$type.1.error.graph,
  size_control_res_3$type.1.error.graph,
  titles = c("γ = 0, (Γ = 1)", 
             "γ = 1, (Γ = 2.718)", 
             "γ = 1.5, (Γ = 4.482)")# Adjust this to match your axis text size
)

# To save the plot:
# ggsave("combined_size_control_plots.png", combined_plot_clean, 
#        width = 13, height = 5, dpi = 300)











##################################################################################
#                                                                                #
#       Here starts the experiments and generate the graphs in the manuscript.   #  
#       We include two presentations. First, the type I error rate at different  #
#       Gamma/gamma's fixing the test statistics, and the margins.               #  
#      or, the type I error rate at a fixing gamma's, fixing the test statistics,#
#      but increase the margins proportionally                                   #
#                                                                                #
##################################################################################



### Treatment margins (60,10,20), outcome margins (15,75), delta = (0,0,1), gamma = 0
size_control_res_1 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(60,10,20), 
                                                                                  outcome.margins = c(15, 75), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 0, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_1$test.stat.hist
size_control_res_1$type.1.error.graph




### Treatment margins (60,10,20), outcome margins (15,75), delta = (0,0,1), gamma = 1
size_control_res_2 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(60,10,20), 
                                                                                  outcome.margins = c(15, 75), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size=18)

size_control_res_2$test.stat.hist
size_control_res_2$type.1.error.graph



### Treatment margins (60,10,20), outcome margins (15,75), delta = (0,0,1), gamma = 1.5
size_control_res_3 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(60,10,20), 
                                                                                  outcome.margins = c(15, 75), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1.5, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_3$test.stat.hist
size_control_res_3$type.1.error.graph






### Treatment margins (30,30,40), outcome margins (15,85), delta = (0,0,1), gamma = 0
size_control_res_4 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(30,30,40), 
                                                                                  outcome.margins = c(15, 85), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 0, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000,base_size = 18)

size_control_res_4$test.stat.hist
size_control_res_4$type.1.error.graph


### Treatment margins (30,30,40), outcome margins (15,85), delta = (0,0,1), gamma = 1.0
size_control_res_5 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(30,30,40), 
                                                                                  outcome.margins = c(15, 85), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1.0, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_5$test.stat.hist
size_control_res_5$type.1.error.graph




### Treatment margins (30,30,40), outcome margins (15,85), delta = (0,0,1), gamma = 1.5
size_control_res_6 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(30,30,40), 
                                                                                  outcome.margins = c(15, 85), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1.5, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_6$test.stat.hist
size_control_res_6$type.1.error.graph





### Treatment margins (60,20,2), outcome margins (10,72), delta = (0,0,1), gamma = 1.0
size_control_res_7 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(60,20,2), 
                                                                                  outcome.margins = c(10, 72), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_7$test.stat.hist
size_control_res_7$type.1.error.graph



### Treatment margins (120,40,4), outcome margins (20,144), delta = (0,0,1), gamma = 1.0
size_control_res_8 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(120,40,4), 
                                                                                  outcome.margins = c(20, 144), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_8$test.stat.hist
size_control_res_8$type.1.error.graph




### Treatment margins (240,80,8), outcome margins (40,288), delta = (0,0,1), gamma = 1.0
size_control_res_9 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(240,80,8), 
                                                                                  outcome.margins = c(40, 288), 
                                                                                  delta = c(0,0,1), 
                                                                                  gamma = 1, 
                                                                                  treatment.scores = c(0,1,2), 
                                                                                  outcome.scores = c(0,1), 
                                                                                  mc.iteration = 1000, base_size = 18)

size_control_res_9$test.stat.hist
size_control_res_9$type.1.error.graph


### Treatment margins (480,160,16), outcome margins (80,576), delta = (0,0,1), gamma = 1.0
size_control_res_10 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(480,160,16), 
                                                                                   outcome.margins = c(80, 576), 
                                                                                   delta = c(0,0,1), 
                                                                                   gamma = 1, 
                                                                                   treatment.scores = c(0,1,2), 
                                                                                   outcome.scores = c(0,1), 
                                                                                   mc.iteration = 1000,base_size = 18)

size_control_res_10$test.stat.hist
size_control_res_10$type.1.error.graph



##### create the combined graph 
# combined_plot_clean <- combine_plots_patchwork(
#  size_control_res_7$type.1.error.graph,
#  size_control_res_9$type.1.error.graph,
#  size_control_res_10$type.1.error.graph,
#  titles = c("(60,20,2); (10,72)", 
#             "(240,80,8); (40,288)", 
#             "(480, 160, 16); (80,576)")# Adjust this to match your axis text size
#) 




### Treatment margins (20,10,10), outcome margins (10,30), delta = (0,0,1), gamma = 1.0
size_control_res_11 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(20,10,10), 
                                                                                   outcome.margins = c(10, 30), 
                                                                                   delta = c(0,0,1), 
                                                                                   gamma = 1, 
                                                                                   treatment.scores = c(0,1,2), 
                                                                                   outcome.scores = c(0,1), 
                                                                                   mc.iteration = 1000,base_size = 18)

size_control_res_11$test.stat.hist
size_control_res_11$type.1.error.graph




### Treatment margins (40,20,20), outcome margins (20,60), delta = (0,0,1), gamma = 1.0
size_control_res_12 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(40,20,20), 
                                                                                   outcome.margins = c(20, 60), 
                                                                                   delta = c(0,0,1), 
                                                                                   gamma = 1, 
                                                                                   treatment.scores = c(0,1,2), 
                                                                                   outcome.scores = c(0,1), 
                                                                                   mc.iteration = 1000,base_size = 18)

size_control_res_12$test.stat.hist
size_control_res_12$type.1.error.graph



### Treatment margins (80,40,40), outcome margins (40,120), delta = (0,0,1), gamma = 1.0
size_control_res_13 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(80,40,40), 
                                                                                   outcome.margins = c(40, 120), 
                                                                                   delta = c(0,0,1), 
                                                                                   gamma = 1, 
                                                                                   treatment.scores = c(0,1,2), 
                                                                                   outcome.scores = c(0,1), 
                                                                                   mc.iteration = 1000,base_size = 18)

size_control_res_13$test.stat.hist
size_control_res_13$type.1.error.graph






### Treatment margins (200,100,100), outcome margins (100,300), delta = (0,0,1), gamma = 1.0
size_control_res_14 = normal.exact.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(200,100,100), 
                                                                                   outcome.margins = c(100, 300), 
                                                                                   delta = c(0,0,1), 
                                                                                   gamma = 1, 
                                                                                   treatment.scores = c(0,1,2), 
                                                                                   outcome.scores = c(0,1), 
                                                                                   mc.iteration = 1000,base_size = 18)

size_control_res_14$test.stat.hist
size_control_res_14$type.1.error.graph


#combined_plot_clean <- combine_plots_patchwork(
#  size_control_res_11$type.1.error.graph,
#  size_control_res_13$type.1.error.graph,
#  size_control_res_14$type.1.error.graph,
#  titles = c("(20,10,10); (10,30)", 
#             "(80,40,40); (40,120)", 
#             "(200, 100, 100); (100,300)")# Adjust this to match your axis text size
#)










################################################################################
#
#
#     When the sample size is large so that exact computation is not feasible,
#     here we present the type I error rate coming from the normal approximation
#     only
#
################################################################################


### A function to sample the data and then compute the type I error rate
### based on the normal approximation 

normal.multiple.treatment.binary.outcome.size.simu.fun = function(
    treatment.margins, 
    outcome.margins, 
    delta, 
    gamma, 
    treatment.scores, 
    outcome.scores, 
    mc.iteration = 10000,
    verbose = FALSE,
    seed = 100,
    base_size = 20
){
  set.seed(seed)
  
  # Check that we have at least 2 treatments
  if(length(treatment.margins) < 2){
    stop("this function requires at least two treatments")
  }
  
  # Check that outcome is binary
  if(length(outcome.margins) != 2){
    stop("this function only discusses binary outcome")
  }
  
  # Check that delta has the same length as treatment.margins
  if(length(delta) != length(treatment.margins)){
    stop("delta should have the same length as treatment.margins")
  }
  
  # Check that treatment.scores has the same length as treatment.margins
  if(length(treatment.scores) != length(treatment.margins)){
    stop("treatment.scores should have the same length as treatment.margins")
  }
  
  # Check that delta is nondecreasing
  if(any(diff(delta) < 0)){
    stop("delta should be nondecreasing")
  }
  
  # Check that treatment.scores is nondecreasing
  if(any(diff(treatment.scores) < 0)){
    stop("treatment.scores should be nondecreasing")
  }
  
  if(sum(treatment.margins) != sum(outcome.margins)){
    stop("what is the correct sample size?")
  }
  
  transform.fun = function(tb){
    scores = outer(treatment.scores, outcome.scores)
    return(sum(tb * scores))
  }
  
  weights = exp(gamma * delta)
  hyper.obs.stat.vec = rep(NA, mc.iteration)
  
  
  for(mc in 1:mc.iteration){
    
    # Sample allocation of outcome 2 across treatments:
    trt.outcome.2 <- as.vector(rMFNCHypergeo(
      nran = 1,
      m = treatment.margins,
      n = outcome.margins[2], # total outcome=2
      odds = weights
    ))
    
    # Outcome 1 by subtraction:
    trt.outcome.1 <- treatment.margins - trt.outcome.2
    
    obs.table <- cbind(trt.outcome.1, trt.outcome.2)
    
    ## compute the test statistic
    hyper.obs.stat.vec[mc] = transform.fun(obs.table)
    
  }
  ## after running the first set of monte carlo, we obtain an estimate of the worst.case.mean and the worst.case.sd
  worst.case.mean = mean(hyper.obs.stat.vec)
  worst.case.sd   =   sd(hyper.obs.stat.vec)
  
  ## now starts another round of monte carlo simulation, but only for pvalue
  normal.p.vec= rep(NA, mc.iteration)
  ## renew the hyper.obs.stat.vec
  hyper.obs.stat.vec=rep(NA,mc.iteration)
  
  
  for(mc in 1:mc.iteration){
    ## use the same sampling steps 
    trt.outcome.2 <- as.vector(rMFNCHypergeo(
      nran = 1,
      m = treatment.margins,
      n = outcome.margins[2], # total outcome=2
      odds = weights
    ))
    
    trt.outcome.1 <- treatment.margins - trt.outcome.2
    
    obs.table <- cbind(trt.outcome.1, trt.outcome.2)
    
    obs.stat  = transform.fun(obs.table)
    ## compute the test statistic
    hyper.obs.stat.vec[mc] = obs.stat 
    
    
    
    ## make sure that both the exact and the normal uses the same worst-case unmeasured confounder 
    normal.p.vec[mc] = norm.score.sen.IxJ(obs.table = obs.table, gamma = gamma, delta = delta, treatment.scores = treatment.scores, outcome.scores = outcome.scores, u_space = t(as.matrix(c(0, outcome.margins[2]))))$max.prob
    
    if(mc %% 10 == 0){
      cat("current mc", mc, "\n")
    }
    
    
  }
  
  ## after all the running, return
  test.stat.hist = hist(hyper.obs.stat.vec)
  
  ## construct the type I error rate graph 
  alpha_seq <- seq(0, 1, length.out = mc.iteration)
  
  df_normal <- data.frame(
    alpha = alpha_seq,
    type1err = sapply(alpha_seq, function(a) mean(normal.p.vec <= a)),
    Method = "Normal"
  )
  
  df_uniform <- data.frame(
    alpha = alpha_seq,
    type1err = alpha_seq,
    Method = "Uniform"
  )
  
  ## report the max of the inflation of type I error rate, which is type1err - alpha_seq
  normal.type.1.err = sapply(alpha_seq, function(a) mean(normal.p.vec <=a))
  # Combine data frames
  plot_df <- rbind(df_normal, df_uniform)
  plot_df$Method <- factor(plot_df$Method, levels = c("Normal", "Uniform"))
  
  # Plotting
  p <- ggplot(plot_df, aes(x = alpha, y = type1err, color = Method, linetype = Method)) +
    geom_line(size = 1.4) +
    scale_color_manual(values = c("#00BFC4", "black")) +
    scale_linetype_manual(values = c("longdash", "dashed")) +
    labs(x = "Nominal Significance Level", y = "Empirical Rejection Rate") +
    theme_minimal(base_size = base_size, base_family = "Times New Roman") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = base_size, family = "Times New Roman"),
      axis.title = element_text(size = base_size, family = "Times New Roman"),
      axis.text = element_text(size = base_size, family = "Times New Roman"),
      legend.key.width = unit(2.5, "cm"),
      legend.position = "bottom"
    )
  return(list(obs.stat.vec = hyper.obs.stat.vec, test.stat.hist = test.stat.hist, type.1.error.graph = p, worst.case.mean = worst.case.mean, worst.case.sd = worst.case.sd, normal.p.vec = normal.p.vec))
  
  
}




## Treatment margins (35,225,240), outcome margins (25,475), delta = (0,0,1), gamma = 0
size_control_res_15 = normal.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(35,225,240), 
                                                                             outcome.margins = c(25, 475), 
                                                                             delta = c(0,0,1), 
                                                                             gamma = 0, 
                                                                             treatment.scores = c(0,1,2), 
                                                                             outcome.scores = c(0,1), 
                                                                             mc.iteration = 1000,base_size = 18)

size_control_res_15$test.stat.hist
size_control_res_15$type.1.error.graph





###  Treatment margins (350,2250,2400), outcome margins (250,4750), delta = (0,0,1), gamma = 0
size_control_res_16 = normal.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(35,225,240)*10, 
                                                                             outcome.margins = c(25, 475)*10, 
                                                                             delta = c(0,0,1), 
                                                                             gamma = 1, 
                                                                             treatment.scores = c(0,1,2), 
                                                                             outcome.scores = c(0,1), 
                                                                             mc.iteration = 1000,base_size = 18)

size_control_res_16$test.stat.hist
size_control_res_16$type.1.error.graph





## Treatment margins (35,225,240), outcome margins (25,475), delta = (0,0,1), gamma = 1.5
size_control_res_17 = normal.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(35,225,240), 
                                                                             outcome.margins = c(25, 475), 
                                                                             delta = c(0,0,1), 
                                                                             gamma = 1.5, 
                                                                             treatment.scores = c(0,1,2), 
                                                                             outcome.scores = c(0,1), 
                                                                             mc.iteration = 1000,base_size = 18)

size_control_res_17$test.stat.hist
size_control_res_17$type.1.error.graph







## Treatment margins (350,2250,2400), outcome margins (250,4750), delta = (0,0,1), gamma = 1.5
size_control_res_18 = normal.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(35,225,240)*10, 
                                                                             outcome.margins = c(25, 475)*10, 
                                                                             delta = c(0,0,1), 
                                                                             gamma = 1.5, 
                                                                             treatment.scores = c(0,1,2), 
                                                                             outcome.scores = c(0,1), 
                                                                             mc.iteration = 1000,base_size = 18)

size_control_res_18$test.stat.hist
size_control_res_18$type.1.error.graph




## Treatment margins (110,195,195), outcome margins (230,270), delta = (0,0,1), gamma = 1.0
size_control_res_19 = normal.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(110,195,195), 
                                                                             outcome.margins = c(230, 270), 
                                                                             delta = c(0,0,1), 
                                                                             gamma = 1, 
                                                                             treatment.scores = c(0,1,2), 
                                                                             outcome.scores = c(0,1), 
                                                                             mc.iteration = 1000,base_size = 18)

size_control_res_19$test.stat.hist
size_control_res_19$type.1.error.graph







## Treatment margins (1100,1950,1950), outcome margins (2300,2700), delta = (0,0,1), gamma = 1.0
size_control_res_20 = normal.multiple.treatment.binary.outcome.size.simu.fun(treatment.margins = c(110,195,195)*10, 
                                                                             outcome.margins = c(230, 270)*10, 
                                                                             delta = c(0,0,1), 
                                                                             gamma = 1, 
                                                                             treatment.scores = c(0,1,2), 
                                                                             outcome.scores = c(0,1), 
                                                                             mc.iteration = 1000,base_size = 18)

size_control_res_20$test.stat.hist
size_control_res_20$type.1.error.graph











