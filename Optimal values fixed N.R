library(tidyverse)
library(directlabels)
library(ggrepel)

#This function here computes the analytical solution to find the minimum for the multiple-objective optimal design
#in a SMART design, under a fixed sample size, using the derived formulae for the computation of the optimal p1, p2
#and p3. The user has to define the following quantities:response rates for the two first-stage treatments (gamma1 
#and gamma2), the weights for the multiple-objective optimal design (w13, w14, w23, w24), such that their sum has 
#to be equal to 1, and the total sample size. The function gives as output the optimal design, 
#relative efficiency of the balanced design, and contour plots of the probabilities p1, p2 and p3 when one 
#at the time is fixed to the optimal value.

optimal <- function(g1 = 0.50, g2 = 0.50, w13 = 0.25, w14 = 0.25, w23 = 0.25, w24 = 0.25, N = 1000){
  time_in <- Sys.time()
  
  #Function value to be minimized. It depends on the constants to be defined by the user and on p1, p2 and p3.
  fun <- function(p1, p2, p3){
    (((g1*p2 + (1-g1))/(N*p1*p2)) + ((g2*p3 + (1-g2))/(N*(1-p1)*p3)))*w13+
      (((g1*p2 + (1-g1))/(N*p1*p2)) + ((g2*(1-p3) + (1-g2))/(N*(1-p1)*(1-p3))))*w14+
      (((g1*(1-p2) + (1-g1))/(N*p1*(1-p2))) + ((g2*p3 + (1-g2))/(N*(1-p1)*p3)))*w23+
      (((g1*(1-p2) + (1-g1))/(N*p1*(1-p2))) + ((g2*(1-p3) + (1-g2))/(N*(1-p1)*(1-p3))))*w24
  }
  
  #Optimal p2, p3 and p1 are computed using their respective formulae
  p2.opt <- sqrt(w13 + w14) / (sqrt(w13 + w14) + sqrt(w23 + w24))
  p3.opt <- sqrt(w13 + w23) / (sqrt(w13 + w23) + sqrt(w14 + w24))
  A <- ((g2*p3.opt + (1-g2))*(w13+w23)*(1-p3.opt) + (g2*(1-p3.opt)+(1-g2))*(w14+w24)*p3.opt)*p2.opt*(1-p2.opt)
  B <- ((g1*p2.opt + (1-g1))*(w13+w14)*(1-p2.opt) + (g1*(1-p2.opt)+(1-g1))*(w23+w24)*p2.opt)*p3.opt*(1-p3.opt)
  p1.opt <- sqrt(B) / (sqrt(A) + sqrt(B))
  
  #The optimal design is computed, using the optimal p1, p2 and p3 derived above.
  opt <- data.frame(g1 = g1, g2 = g2, p1.opt = round(p1.opt,4), p2.opt = round(p2.opt,4), p3.opt = round(p3.opt,4), 
                    fun = fun(p1.opt, p2.opt, p3.opt),N = N, w13 = w13, w14 = w14, w23 = w23, w24 = w24)
  
  #Relative efficiency and effective sample size of the balanced design are computed
  r.eff <- opt$fun / fun(0.50, 0.50, 0.50)
  N.eff <- N * (1/r.eff)
  
  #Plot relative efficiency
  
  #A dataframe including possible combinations of p1, p2 and p3 is created. Then relative efficiency with respect
  #to the optimal design is computed
  re_data <- data.frame(expand.grid(p1 = seq(0.01, 0.99, 0.03), p2 = seq(0.01, 0.99, 0.03), 
                                    p3 = seq(0.01, 0.99, 0.03)), fun = NA, r.eff = NA)
  
  ff <- opt %>% mutate(r.eff = 1) %>% select(p1.opt, p2.opt, p3.opt, fun, r.eff)
  colnames(ff) <- colnames(re_data)
  
  re_data <- re_data %>% rowwise() %>% mutate(fun = fun(p1, p2, p3), r.eff = round((ff$fun / fun),4)) %>% 
    arrange(fun)

  
  #Contour plots
  #p1 fixed to optimal. p2 and p3 vary across their domain
  cont_data1 <- data.frame(p1 = opt$p1.opt,
                           expand.grid(p2 = seq(0.01, 0.99, 0.01), p3 = seq(0.01, 0.99, 0.01)),
                           fun = NA, r.eff = NA)
  
  cont_data1 <- cont_data1 %>% rowwise() %>% mutate(fun = fun(p1, p2, p3)) %>% 
    arrange(fun)
  cont_data1 <- cont_data1 %>% rowwise() %>% mutate(r.eff = (ff$fun / fun)) 
  
  cont_plot1 <- cont_data1 %>% ggplot(aes(x = p2, y = p3, z = r.eff))+ stat_contour(breaks = c(0.20,0.40,
                                                                                               0.60,0.70,0.80,
                                                                                               0.90,0.95)) + 
    geom_point(data = ff, aes(x = p2, y = p3, colour = "firebrick")) + 
    geom_label_repel(data = ff, aes(x = p2, y = p3, label = "Optimal Design"), vjust = 2.4,
                     segment.size = 0.2) +
    geom_dl(aes(label=..level.., colour = "firebrick"), method=list("last.bumpup",hjust = 1.1, cex = 0.8), 
            stat="contour",breaks = c(0.20,0.40,
                                      0.60,0.70,0.80,
                                      0.90,0.95)) +
    theme_minimal() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank() ,axis.line.x = element_line(colour = "black"),
                            axis.line.y = element_line(colour = "black"), legend.position = "none") +
    ggtitle(expression(paste("Contour graph for ", p[2], " and ", p[3], " - ", p[1], " fixed to ",p[1]^"*")),
            subtitle = paste("Optimal design:","p1* =", paste0(ff$p1,","), 
                             "p2* =", paste0(ff$p2, ","), "p3* =", paste0(ff$p3))) +
    labs(x = expression(p[2]), y = expression(p[3])) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_x_continuous(breaks = seq(0,1,0.1)) +
    scale_y_continuous(breaks = seq(0,1,0.1))
  
#p2 fixed to optimal. p1 and p3 vary across their domain
  cont_data2 <- data.frame(p2 = opt$p2.opt,
                           expand.grid(p1 = seq(0, 1, 0.01), p3 = seq(0, 1, 0.01)),
                           fun = NA, r.eff = NA)
  
  cont_data2 <- cont_data2 %>% rowwise() %>% mutate(fun = fun(p1, p2, p3)) %>% 
    arrange(fun)
  cont_data2 <- cont_data2 %>% rowwise() %>% mutate(r.eff = (ff$fun / fun))
  
  cont_plot2 <- cont_data2 %>% ggplot(aes(x = p1, y = p3, z = r.eff)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                0.60,0.70,0.80,
                                                                                                0.90,0.95)) +
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
  
  #p3 fixed to optimal. p1 and p2 vary across their domain
  cont_data3 <- data.frame(p3 = opt$p3.opt,
                           expand.grid(p1 = seq(0, 1, 0.01), p2 = seq(0, 1, 0.01)),
                          fun = NA, r.eff = NA)
  
  cont_data3 <- cont_data3 %>% rowwise() %>% mutate(fun = fun(p1, p2, p3)) %>% 
    arrange(fun)
  cont_data3 <- cont_data3 %>% rowwise() %>% mutate(r.eff = (ff$fun / fun))
  
  cont_plot3 <- cont_data3 %>% ggplot(aes(x = p1, y = p2, z = r.eff)) + stat_contour(breaks = c(0.20,0.40,
                                                                                                0.60,0.70,0.80,
                                                                                                0.90,0.95)) +
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
    
  
  colnames(opt) <- c(paste0(intToUtf8(947),c(1:2)), "p1.opt", "p2.opt", "p3.opt",
                       paste(intToUtf8(966), "value"), "N", paste0(intToUtf8(955),c("13","14","23","24")))
  
  time_fin <- Sys.time() - time_in
  list(opt = opt, r.eff = r.eff, N.eff = N.eff, time = time_fin,
       cont_plot1 = cont_plot1, cont_plot2 = cont_plot2, cont_plot3 = cont_plot3)
}
