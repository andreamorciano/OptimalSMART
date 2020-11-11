library(tidyverse)
library(ggrepel)
library(directlabels)

#This function here computes the numerical solution to find the minimum for the multiple-objective optimal design
#in a SMART design, when considering cost of treatments into account. The solution is given by using an iterative
#algorithm that evaluates the function of interest under different combinations of the variables p1, p2 and p3,
#and chooses the one that provides the lowest function value. The user has to define the following quantities:
#response rates for the two first-stage treatments (gamma1 and gamma2), cost of treatments (Cp and Cn), the weights
#for the multiple-objective optimal design (w13, w14, w23, w24), that their sum has to be equal to 1, and the total
#budget. The function gives as output the optimal design, relative efficiency of the balanced design, and contour plots
#of the probabilities p1, p2 and p3 when one at the time is fixed to the optimal value.

analytic.prop <- function(gamma1 = 0.5, gamma2 = 0.5, Ca = 10, Cb = 10, Cc = 10, Cd = 10,
                          Ce = 10, Cf = 10, Cg = 10, Ch = 10,w13 = 0.25, w23 = 0.25, w14 = 0.25, 
                          w24 = 0.25, C = 100000){
  time_in <- Sys.time()
  g1 <- gamma1
  g2 <- gamma2
  
  #function stops if the sum of weights is different from 1
  wgts <- w13 + w14 + w23 + w24
  stopifnot(near(wgts, 1))
  
  
  #Function to compute sample size, given p1, p2, p3 and the constants.
  N.opt <- function(x1, x2, x3){
    N <- C / (x1*(Ca+Cc*g1+(1-g1)*x2*Cd+(1-g1)*(1-x2)*Ce) + (1-x1)*(Cb+Cf*g2+(1-g2)*x3*Cg+(1-g2)*(1-x3)*Ch))
    N <- round(N, digits = 9) 
  }
  
  ##Function that computes total cost
  Cost <- function(N, x1, x2, x3){
    C <- N*(x1*(Ca+Cc*g1+(1-g1)*x2*Cd+(1-g1)*(1-x2)*Ce) + (1-x1)*(Cb+Cf*g2+(1-g2)*x3*Cg+(1-g2)*(1-x3)*Ch))
  }
  
  #Function value, given p1, p2, p3, N and the constants.
  fun <- function(x1, x2, x3, N){
    ((g1*x2+(1-g1))/(N*x1*x2) + (g2*x3+(1-g2))/(N*(1-x1)*x3))*w13 +
      ((g1*x2+(1-g1))/(N*x1*x2) + (g2*(1-x3)+(1-g2))/(N*(1-x1)*(1-x3)))*w14 +
      ((g1*(1-x2)+(1-g1))/(N*x1*(1-x2)) + (g2*x3+(1-g2))/(N*(1-x1)*x3))*w23 +
      ((g1*(1-x2)+(1-g1))/(N*x1*(1-x2)) + (g2*(1-x3)+(1-g2))/(N*(1-x1)*(1-x3)))*w24
  }
  
  #The temporary minimum of the function for the evaluated combinations of p1, p2 and p3 is computed
  iterative <- function(p1, p2, p3){
    #A dataframe with possible combinations of p1, p2 and p3 is created and the function value, the sample size
    #and the cost are computed. Then the combination that provides the minimum value is chosen as the local temporary minimum.
    values <- data.frame(expand.grid(p1 = p1, p2 = p2, p3 = p3), N = NA, fun = NA, Cost = NA)
    values <- values %>% rowwise() %>% mutate(N = N.opt(p1, p2, p3), fun = fun(p1, p2, p3, N),Cost = Cost(N, p1, p2, p3))
    intermed <- values %>% arrange(fun) %>% filter(fun > 0 & p1>0 & p1<1 & p2<1 & p2>0 & p3<1 & p3>0) %>% slice(which.min(fun))
    return(intermed)
  }
  
  #p1, p2 and p3 are initialized
  p1 <- p2 <- p3 <- seq(0.01, 0.99, by = 0.02)
  
  #First iteration is computed
  iterations <- data.frame(1, iterative(p1,p2,p3))
  colnames(iterations) <- c("Iteration #", "p1", "p2", "p3", "N", "fun", "Cost")
  
  i <- 2
  
  #The step size and the size of the interval change according to the iteration number. For iteration #2 the algorithm
  #is distint than for the other iterations.
  while (i < 100) {
    intermed <- iterations[i-1, -1]
    
    if(i != 2){
      k <- 3.5 * 10^(-(i-1))
      
      p1 <- seq((intermed$p1 - k), (intermed$p1 + k), by = 2 * 10^(-i))
      p2 <- seq((intermed$p2 - k), (intermed$p2 + k), by = 2 * 10^(-i))
      p3 <- seq((intermed$p3 - k), (intermed$p3 + k), by = 2 * 10^(-i))
    } else {
      k <- 1.8 * 10^(-(i-1))
      
      p1 <- seq((intermed$p1 - k), (intermed$p1 + k), by = 0.3* 10^(-i+1))
      p2 <- seq((intermed$p2 - k), (intermed$p2 + k), by = 0.3 * 10^(-i+1))
      p3 <- seq((intermed$p3 - k), (intermed$p3 + k), by = 0.3 * 10^(-i+1))
    }
    
    temp <- data.frame(i, iterative(p1,p2,p3))
    colnames(temp) <- c("Iteration #", "p1", "p2", "p3", "N", "fun", "Cost")
    
    iterations <- rbind(iterations, temp)
    
    #Convergence criteria that is met if the results found in two subsequent iterations are the same up to a certain
    #degree of precision
    
    i <- ifelse(near(iterations[i-1,2], iterations[i,2], tol = .Machine$double.eps^0.3)
                & near(iterations[i-1,3], iterations[i,3], tol = .Machine$double.eps^0.3) & 
                  near(iterations[i-1,4], iterations[i,4], tol = .Machine$double.eps^0.3) & i > 5, 100, i + 1)
  }
  
  #Optimal design is found after the convergence criterion is met
  final <- iterations %>% transmute(gamma1 = g1, gamma2 = g2, Ca = Ca, Cb = Cb, Cc = Cc, Cd = Cd, Ce = Ce, Cf = Cf,
                                    Cg = Cg, Ch = Ch,
                                    p1.opt = round(p1, 4), p2.opt = round(p2, 4), p3.opt = round(p3, 4),
                                    N.opt = round(N,4), fun = round(fun, 8), Cost = round(Cost),
                                    w13 = w13, w14 = w14, w23 = w23, w24 = w24) %>% slice(which.max(iterations$`Iteration #`))
  
  #Relative efficiency is computed as the ratio of the function value of the balanced design and that of the optimal
  #design found above. Cost efficiency is then computed
  
  #Balanced design
  N_bal <- N.opt(x1 = 0.5, x2 = 0.5, x3 = 0.5)
  fun_bal <- fun(x1 = 0.5, x2 = 0.5, x3 = 0.5, N = N_bal)
  
  #Relative efficiency and cost efficiency
  r.eff <- final$fun / fun_bal
  cost.eff <- final$Cost * (1/r.eff)
  
  #Plot relative efficiency
  #A dataframe with possible combinations of p1, p2 and p3 is created and relative efficiency with respect to the
  #optimal design is computed
  
  re_data <- data.frame(expand.grid(p1 = seq(0.01, 0.99, 0.03), p2 = seq(0.01, 0.99, 0.03), 
                                    p3 = seq(0.01, 0.99, 0.03)), N = NA, fun = NA, r.eff = NA, eff_cost = NA)
  ff <- final %>% mutate(r.eff = 1) %>% select(p1.opt, p2.opt, p3.opt, N.opt, fun, r.eff, Cost)
  colnames(ff) <- colnames(re_data)
  
  re_data <- re_data %>% rowwise() %>% mutate(N = N.opt(p1, p2, p3), fun = fun(p1, p2, p3, N), r.eff = round((ff$fun / fun),4), eff_cost = round(C * (1/r.eff),4)) %>% 
    arrange(fun)

  #Column names for the object containing the optimal design are set
  colnames(final) <- c(paste0(intToUtf8(947),c(1:2)), "Ca", "Cb", "Cc", "Cd", "Ce", "Cf", "Cg", "Ch",
                       "p1.opt", "p2.opt", "p3.opt", "N.opt",
                       paste(intToUtf8(966), "value"), "Cost", paste0(intToUtf8(955),c("13","14","23","24")))
  
  #Contour plots
  #p1 fixed to its optimal value. Contour plots are for p2 and p3.
  cont_data1 <- data.frame(p1 = final$p1.opt,
                           expand.grid(p2 = seq(0.01, 0.99, 0.01), p3 = seq(0.01, 0.99, 0.01)),
                           N = NA, fun = NA, r.eff = NA, eff_cost = NA)
  
  cont_data1 <- cont_data1 %>% rowwise() %>% mutate(N = N.opt(p1, p2, p3), fun = fun(p1, p2, p3, N)) %>% 
    arrange(fun)
  cont_data1 <- cont_data1 %>% rowwise() %>% mutate(r.eff = (ff$fun / fun) , eff_cost = C * (1/r.eff))

  cont_plot1 <- cont_data1 %>% ggplot(aes(x = p2, y = p3, z = r.eff)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                0.60,0.70,0.80,
                                                                                                0.90,0.95,1)) +
    geom_point(data = ff, aes(x = p2, y = p3, colour = "firebrick")) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) + 
    geom_label_repel(data = ff, aes(x = p2, y = p3, label = "Optimal Design"), vjust = 2.4,
                     segment.size = 0.2) + 
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[2], " and ", p[3], " - ", p[1], " fixed to ",p[1]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(ff$p1,","), 
                             "p2* =", paste0(ff$p2, ","), "p3* =", paste0(ff$p3))) +
    labs(x = expression(p[2]), y = expression(p[3])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))
  
  #p2 fixed to its optimal value. contour plots are for p1 and p3 varying across their domain.
  cont_data2 <- data.frame(p2 = final$p2.opt,
                           expand.grid(p1 = seq(0.01, 0.99, 0.01), p3 = seq(0.01, 0.99, 0.01)),
                           N = NA, fun = NA, r.eff = NA, eff_cost = NA)
  
  cont_data2 <- cont_data2 %>% rowwise() %>% mutate(N = N.opt(p1, p2, p3), fun = fun(p1, p2, p3, N)) %>% 
    arrange(fun)
  cont_data2 <- cont_data2 %>% rowwise() %>% mutate(r.eff = (ff$fun / fun) , eff_cost = C * (1/r.eff)) 
  
  cont_plot2 <- cont_data2 %>% ggplot(aes(x = p1, y = p3, z = r.eff)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                0.60,0.70,0.80,
                                                                                                0.90,0.95,1)) +
    geom_point(data = ff, aes(x = p1, y = p3, colour = "firebrick")) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) + 
    geom_label_repel(data = ff, aes(x = p1, y = p3, label = "Optimal Design"), vjust = 2.4,
                     segment.size = 0.2) + 
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[1], " and ", p[3], " - ", p[2], " fixed to ",p[2]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(ff$p1,","), 
                             "p2* =", paste0(ff$p2, ","), "p3* =", paste0(ff$p3))) +
    labs(x = expression(p[1]), y = expression(p[3])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))
  
  #p3 fixed to its optimal value. p1 and p2 are free to vary across their domain.
  cont_data3 <- data.frame(p3 = final$p3.opt,
                           expand.grid(p1 = seq(0.01, 0.99, 0.01), p2 = seq(0.01, 0.99, 0.01)),
                           N = NA, fun = NA, r.eff = NA, eff_cost = NA)
  
  cont_data3 <- cont_data3 %>% rowwise() %>% mutate(N = N.opt(p1, p2, p3), fun = fun(p1, p2, p3, N)) %>% 
    arrange(fun)
  cont_data3 <- cont_data3 %>% rowwise() %>% mutate(r.eff = (ff$fun / fun) , eff_cost = C * (1/r.eff))
  
  cont_plot3 <- cont_data3 %>% ggplot(aes(x = p1, y = p2, z = r.eff)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                0.60,0.70,0.80,
                                                                                                0.90,0.95,1)) +
    geom_point(data = ff, aes(x = p1, y = p2, colour = "firebrick")) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) + 
    geom_label_repel(data = ff, aes(x = p1, y = p2, label = "Optimal Design"), vjust = 2.4,
                     segment.size = 0.2) + 
    theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[1], " and ", p[2], " - ", p[3], " fixed to ",p[3]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(ff$p1,","), 
                              "p2* =", paste0(ff$p2, ","), "p3* =", paste0(ff$p3))) +
    labs(x = expression(p[1]), y = expression(p[2])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))

  #Time elapsed
  time_el <- Sys.time() - time_in
  
  list(final = final, iterations = iterations, time = time_el, r.eff = r.eff, cost.eff = cost.eff,
       cont_plot1 = cont_plot1, cont_plot2 = cont_plot2, cont_plot3 = cont_plot3)
} 
