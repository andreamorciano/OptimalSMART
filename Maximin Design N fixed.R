library(tidyverse)
library(ggrepel)
library(directlabels)

#This function here computes the maximin design (MMD) as a solution for the misspecification of the response rates
#in a SMART design, under a fixed sample size.  The MMD  is chosen as the design in the design space that maximizes
#the minimum RE over the parameter space, which is the design that performs better in the worst-case scenario.
#The user has to define the following quantities:response rates for the two first-stage treatments (gamma1 
#and gamma2), the weights for the multiple-objective optimal design (w13, w14, w23, w24), such that their sum has 
#to be equal to 1, the total sample size and the stepsize when varying p1, p2 and p3 across their domain.
#An interval for the response rates defines the parameter space for the response rates. This is obtained as an 
#interval with size 0.10 with the point estimate provided by the user being the center of such interval.
#The function gives as output the maximin design, minimum relative efficiency of the balanced design, and 
#contour plots of the probabilities p1, p2 and p3 when one at the time is fixed to its value in the MMD.

maximin_design_nfixed <- function(gamma1, gamma2, w13 = 0.25, w14 = 0.25, w23 = 0.25, w24 = 0.25, N = 1000, stepsize = 0.05){
  time_in <- Sys.time() 
  
  #A priori interval for the response rates is created
  gamma1 <- seq((gamma1 -0.05), (gamma1 + 0.05), 0.01)
  gamma2 <- seq((gamma2 - 0.05), (gamma2 + 0.05), 0.01)
  #The function stops if the weights don't sum up to 1
  wgts <- w13 + w14 + w23 + w24
  stopifnot(near(wgts, 1))
  
  #p1, p2 and p3 are initialized.
  p1 <- p2 <- p3 <- seq(0.01, 0.99, stepsize)

  #Function value, given p1, p2, p3, N and the constants.
  fun <- function(p1, p2, p3){
    (((g1*p2 + (1-g1))/(N*p1*p2)) + ((g2*p3 + (1-g2))/(N*(1-p1)*p3)))*w13+
      (((g1*p2 + (1-g1))/(N*p1*p2)) + ((g2*(1-p3) + (1-g2))/(N*(1-p1)*(1-p3))))*w14+
      (((g1*(1-p2) + (1-g1))/(N*p1*(1-p2))) + ((g2*p3 + (1-g2))/(N*(1-p1)*p3)))*w23+
      (((g1*(1-p2) + (1-g1))/(N*p1*(1-p2))) + ((g2*(1-p3) + (1-g2))/(N*(1-p1)*(1-p3))))*w24
  }
  
  #A dataframe is created, containing all designs in the design space for each combination of the response rates
  #in the parameter space.
  values_mmd <- data.frame(expand.grid(gamma1 = gamma1, gamma2 = gamma2, p1 = p1, p2 = p2, p3 = p3), N = N,
                           fun = NA, RE = NA)
  
  #Design and parameter spaces are created.
  par_space <- expand.grid(gamma1 = gamma1, gamma2 = gamma2)
  des_space <- expand.grid(p1 = p1, p2 = p2, p3 = p3)
  
  #Temporary dataframe created to save the values obtained in the following for loop
  temp <- as.data.frame(matrix(NA,ncol = 8))
  colnames(temp) <- colnames(values_mmd)
  
  #For each combination of the response rates in the parameter space, the locally optimal design is computed. Which is
  #the design in the design space that provides the minimum value for the function of interest.
  for (i in 1:nrow(par_space)) {
    g1 <- par_space[i,1]
    g2 <- par_space[i,2]
    lod <- values_mmd %>% filter(gamma1 == g1 & gamma2 == g2)
    lod <- lod %>% rowwise() %>% mutate(fun = fun(p1, p2, p3)) %>% arrange(fun)
    lod_val <- lod  %>% slice(which.min(fun)) %>% select(fun)
    lod_val <- as.numeric(lod_val)
    lod <- lod %>% rowwise() %>% mutate(RE = lod_val / fun)
    
    temp <- rbind(temp, lod)
  }
  
  values_mmd <- temp %>% filter(!is.na(gamma1))
  
  #Temporary dataframe used to save the results obtained in the following for loop
  temp2 <- as.data.frame(matrix(NA,ncol = 8))
  colnames(temp2) <- c(colnames(values_mmd[1:7]),"minRE")
  
  #For each design in the design space, the minimum relative efficiency is computed, which is the lowest 
  #efficiency that the design achieves all over the parameter space. Then the design having the highest minimum RE
  #is chosen as the MMD.
  for (k in 1:nrow(des_space)) {
    p11 <- des_space[k,1]
    p21 <- des_space[k,2]
    p31 <- des_space[k,3]
    mmd <- values_mmd %>% filter(p1 == p11 & p2 == p21 & p3 == p31) %>% slice(which.min(RE))
    mmd <- mmd %>% mutate(minRE = round(RE,5)) %>% select(-RE)
    
    temp2 <- rbind(temp2,mmd)
  }
  
  g1_ex <- paste0("[",paste(min(gamma1),max(gamma1), sep = ", "),"]")
  g2_ex <- paste0("[",paste(min(gamma2),max(gamma2), sep = ", "),"]")
  
  #MMD is saved in a separate object.
  mmd <- temp2 %>% filter(!is.na(gamma1)) %>% arrange(desc(minRE))
  mmd <- mmd %>% dplyr::transmute(p1_mmd = p1, p2_mmd = p2, p3_mmd = p3, N = N, fun = fun, minRE = minRE)
  mmd_value <- mmd %>% slice(which.max(minRE)) %>% dplyr::transmute(gamma1 = g1_ex, gamma2 = g2_ex, p1_mmd = p1_mmd,
                                                                    p2_mmd = p2_mmd, p3_mmd = p3_mmd, N = N, fun = fun, MMV = minRE)
  
  #Relative efficiency and effective sample size of the balanced design
  r.eff <- mmd %>% filter(p1_mmd == 0.50 & p2_mmd == 0.50 & p3_mmd == 0.50) %>% select(minRE)
  N.eff <- mmd_value$N * (1/r.eff)
  
  time_fin = Sys.time() - time_in
  
  mmd_value <- mmd_value %>% dplyr::transmute(gamma1 = gamma1, gamma2 = gamma2, p1_mmd = p1_mmd, 
                                              p2_mmd = p2_mmd, p3_mmd = p3_mmd, N = N, fun = fun, MMV = MMV,
                                              w13 = w13, w14 = w14, w23 = w23, w24 = w24)
  
  colnames(mmd_value) <- c(paste0(intToUtf8(947),c(1:2)," Range"), "p1.mmd", "p2.mmd", "p3.mmd", 
                           "N", paste(intToUtf8(966), "value"), "MMV", 
                           paste0(intToUtf8(955),c("13","14","23","24")))
  
  #Contour plots for p1, p2 and p3 when fixing one at the time to its value in the MMD design.
  mmd_value2 <- mmd_value %>% mutate(minRE = MMV) %>% select(-MMV)
  
  #Contour plot for p2 and p3. p1 fixed
  cont_data1 <- mmd %>% filter(p1_mmd == mmd_value$p1.mmd)
  
  cont_plot1 <- cont_data1 %>% ggplot(aes(x = p2_mmd, y = p3_mmd, z = minRE)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                        0.60,0.70,0.80,
                                                                                                        0.90,0.95)) +
    geom_point(data = mmd_value2, aes(x = p2.mmd, y = p3.mmd, colour = "firebrick")) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) +
    geom_label_repel(data = mmd_value2, aes(x = p2.mmd, y = p3.mmd, label = "Maximin Design"), vjust = 2.4,
                     segment.size = 0.2) + 
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[2], " and ", p[3], " - ", p[1], " fixed to ",p[1]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(mmd_value2$p1.mmd,","), 
                             "p2* =", paste0(mmd_value2$p2.mmd, ","), "p3* =", paste0(mmd_value2$p3.mmd))) +
    labs(x = expression(p[2]), y = expression(p[3])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))
  
  #Contour plot for p1 and p3. p2 fixed
  cont_data2 <- mmd %>% filter(p2_mmd == mmd_value$p2.mmd)
  
  cont_plot2 <- cont_data2 %>% ggplot(aes(x = p1_mmd, y = p3_mmd, z = minRE)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                        0.60,0.70,0.80,
                                                                                                        0.90,0.95)) +
    geom_point(data = mmd_value2, aes(x = p1.mmd, y = p3.mmd, colour = "firebrick")) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) + 
    geom_label_repel(data = mmd_value2, aes(x = p1.mmd, y = p3.mmd, label = "Maximin Design"), vjust = 2.4,
                     segment.size = 0.2) + 
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[1], " and ", p[3], " - ", p[2], " fixed to ",p[2]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(mmd_value2$p1.mmd,","), 
                             "p2* =", paste0(mmd_value2$p2.mmd, ","), "p3* =", paste0(mmd_value2$p3.mmd))) +
    labs(x = expression(p[1]), y = expression(p[3])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))
  
  #Contour plot for p1 and p2. p3 fixed
  cont_data3 <- mmd %>% filter(p3_mmd == mmd_value$p3.mmd)
  
  cont_plot3 <- cont_data3 %>% ggplot(aes(x = p1_mmd, y = p2_mmd, z = minRE)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                        0.60,0.70,0.80,
                                                                                                        0.90,0.95)) +
    geom_point(data = mmd_value2, aes(x = p1.mmd, y = p2.mmd, colour = "firebrick")) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) +
    geom_label_repel(data = mmd_value2, aes(x = p1.mmd, y = p2.mmd, label = "Maximin Design"), vjust = 2.4,
                     segment.size = 0.2) + 
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[1], " and ", p[2], " - ", p[3], " fixed to ",p[3]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(mmd_value2$p1.mmd,","), 
                             "p2* =", paste0(mmd_value2$p2.mmd, ","), "p3* =", paste0(mmd_value2$p3.mmd))) +
    labs(x = expression(p[1]), y = expression(p[2])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))
  
  list(values_mmd = values_mmd, mmd_value = mmd_value, mmd = mmd, time = time_fin, r.eff = r.eff,
       N.eff = N.eff, cont_plot1 = cont_plot1, cont_plot2 = cont_plot2, cont_plot3 = cont_plot3)
}
