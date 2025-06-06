# Statistics ------
CI_ME_95 <- function(X,retCI = F){
  n <- length(X)
  xbar <- mean(X)
  sX <- sd(X)
  ME <- qt(0.975,df = n-1)*sX/sqrt(n)
  
  if(retCI){
    low_int <- xbar - ME
    high_int <- xbar + ME
    CI_res <- c(low_int,high_int)
    return(CI_res)
  }else{
    return(ME)
  }
  
}
calc_stats <- function(x, probs = c(0.25, 0.5, 0.75, 0.95)) {
  tibble(
    val = c(mean(x,na.rm=T),sd(x,na.rm=T),CI_ME_95(x),median(x,na.rm=T),max(x,na.rm=T),quantile(x, probs, na.rm = TRUE)),
    statistic = c('mean','stdev','CI95_ME','median','maximum',paste(probs,'_quant',sep=''))
  )
}

# Mixing Ratio to Mass Density ------

px_Cx_df <- function(df,R = 8.314472){
  
  ### Testing ###
  # df <- hr_data
  # R <- 8.314472 # m3 Pa (K mol)^-1
  ### Testing ###
  
  ## Interpolate Across Gaps
  df <- df %>% mutate_at(vars(T_C,P_Pa),~imputeTS::na_interpolation(.))
  ## Convert to mass density
  df_px <- df %>% 
    mutate_at(vars(Concentration),~ifelse(!Species %in% c('BC','Org'),
                                          ((.*P_Pa*MW)/(R*(T_C+273.15)))*10^6,.))
  
  return(df_px)
}

# Image to ggplot2 object -------
img2gg <- function(img_path){
  
  # Read in the image
  img <- image_read(img_path)
  
  # Convert the image to a ggplot object
  ggimg <- ggplot() +
    annotation_custom(
      rasterGrob(img),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    theme_void()
  
  return(ggimg)
}
# Insert NAs for time gaps -------------
ts_gap <- function(df,time_col = 'Time_AKST',gap_val = NA,gap_by = '1 hour'){
  
  ### Testing ###
  # df = hr_data
  # time_col = 'Time_AKST'
  # gap_val = NA
  # gap_by = '1 hour'
  ### Testing ###
  
  time_df <- data.frame(time_col = seq.POSIXt(min(df[[time_col]]),max(df[[time_col]]),gap_by))
  time_df <- setNames(time_df,time_col)
  
  df_gap <- merge(time_df,df,by=time_col,all = T)
  
  return(df_gap)
}
# Figure 1: Time Series --------

gen_figure_1 <- function(df,CTC_df){
  
  ### Testing ###
  # df <- hr_data
  # CTC_df <- CTC_data
  ### Testing ###
  
  panel_1 <- c('CTC_SO2','CTC_temp_3m')
  panel_2 <- c('CO','NOx','sum_aromatics','sum_furans')
  panel_3 <- c('PM1','SO4')
  
  p1col <- c('#CD2626', '#0000CD')
  p2col <- c('#8B1C62',"#FF7F24",'#4F94CD', '#FFC125')
  p3col <- c('#556B2F', '#EEB422')
  
  df <- df %>% mutate('sum_furans' = (Furan+Methylfuran+Furfural+Methylfurfural)*100) %>%
    mutate('sum_aromatics' = Benzene+Toluene+C8_Aromatics+C9_Aromatics+Naphthalene+C10_Aromatics) %>% 
    mutate_at(vars(CO),~.*100) %>%
    select(Time_AKST,all_of(c(panel_2,panel_3))) %>%
    ts_gap() %>%
    pivot_longer(!Time_AKST,names_to = 'Species',values_to = 'val')
  
  CTC_df <- CTC_df %>% filter(between(Time_AKST,df$Time_AKST[1],df$Time_AKST[nrow(df)])) %>% 
    select(Time_AKST,all_of(panel_1)) %>%
    pivot_longer(!Time_AKST,names_to = 'Species',values_to = 'val')
  
  cold_event <- c(as.POSIXct('2022-01-29 00:00:00',tz = 'America/Sitka'),
                  as.POSIXct('2022-02-03 12:00:00',tz = 'America/Sitka'))
  
  warm_event <- c(as.POSIXct('2022-02-23 00:00:00',tz = 'America/Sitka'),
                  as.POSIXct('2022-02-25 13:00:00',tz = 'America/Sitka'))
  
  plt_1 <- ggplot(CTC_df,aes(x=Time_AKST,y=val,color = Species))+
    annotate("rect", xmin=cold_event[1], xmax=cold_event[2], ymin=-40, ymax=55, alpha=0.3, fill="#1E90FF",color = '#1E90FF')+
    annotate("rect", xmin=warm_event[1], xmax=warm_event[2], ymin=-40, ymax=55, alpha=0.2, fill='orange',color = 'orange')+
    annotate("text",label = 'Cold Event',x=cold_event[1]+days(2)+hours(16),y=50,color = 'black',fontface='bold')+
    annotate("text",label = 'Warm\nEvent',x=warm_event[1]+days(1)+hours(6),y=45,color = 'black',fontface='bold')+
    geom_hline(yintercept = 0,linewidth = 1,color = 'black')+
    geom_line(linewidth=1)+
    scale_y_continuous(name = expression(bold(SO[2]*' (ppbv)')),
                       breaks = seq(-40,55,10),minor_breaks = seq(-40,55,2),limits = c(-40,55),expand = c(0,0),
                       sec.axis = sec_axis(transform = ~.,breaks = seq(-40,55,10),
                                           name=expression(bold('Temperature ('*degree*'C)'))))+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.15,0.95),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p1col,labels = c(expression(bold(SO[2])),expression(bold('Temperature'))))+
    guides(color = guide_legend(nrow = 1),
     # y = guide_axis(minor.ticks = TRUE),
      x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  plt_2 <- ggplot(df %>% filter(Species %in% panel_2),aes(x=Time_AKST,y=val,color = Species))+
    annotate("rect", xmin=cold_event[1], xmax=cold_event[2], ymin=-Inf, ymax=Inf, alpha=0.3, fill='#1E90FF',color = '#1E90FF')+
    annotate("rect", xmin=warm_event[1], xmax=warm_event[2], ymin=-Inf, ymax=Inf, alpha=0.2, fill='orange',color = 'orange')+
    geom_line(linewidth=1)+
    scale_y_continuous(name = expression(bold(NO[x]*' and '*Sigma*'Aromatics (ppbv)')),
                       breaks = seq(0,300,100),minor_breaks = seq(0,300,25),limits = c(0,300),expand = c(0,0),
                       sec.axis = sec_axis(transform = ~./100,breaks = seq(0,3,0.5),
                                           name=expression(bold('CO (ppmv) and '*Sigma*'Furanoids (ppbv)'))))+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.15,0.9),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p2col,
                       labels = c(expression(bold('CO')),expression(bold(NO[x])),
                                  expression(bold(Sigma*'Aromatics')),
                                  expression(bold(Sigma*'Furanoids'))))+
    guides(color = guide_legend(nrow = 2),
           #y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  df <- df %>% mutate_at(vars(val),~ifelse(Species == 'SO4',./0.40,.))
  
  plt_3 <- ggplot(df %>% filter(Species %in% panel_3),aes(x=Time_AKST,y=val,color = Species))+
    annotate("rect", xmin=cold_event[1], xmax=cold_event[2], ymin=-Inf, ymax=Inf, alpha=0.3, fill='#1E90FF',color = '#1E90FF')+
    annotate("rect", xmin=warm_event[1], xmax=warm_event[2], ymin=0, ymax=100, alpha=0.2, fill='orange',color = 'orange')+
    geom_line(linewidth=1)+
    scale_y_continuous(name = expression(bold(PM[1]*' ('*mu*'g '*m^-3*')')),
                       breaks = seq(0,100,20),
                       minor_breaks = seq(0, 100, by = 5),
                       limits = c(0,100),expand = c(0,0),
                       sec.axis = sec_axis(name = expression(bold({SO[4]^{'2-'}}*' ('*mu*'g '*m^-3*')')),
                                           transform = ~.*0.40,breaks = seq(0,40,5)))+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    xlab('Local Time (AKST)')+
    theme(legend.position = 'inside',legend.position.inside = c(0.13,0.95),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p3col,labels = c(expression(bold(PM[1])),expression(bold({SO[4]^{'2-'}}))))+
      guides(color = guide_legend(nrow = 1),
             #y = guide_axis(minor.ticks = TRUE),
             x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  fin_plot <- (plt_1/plt_2/plt_3 & theme(plot.tag.position  = c(0.075, .96)))+
    plot_annotation(tag_levels = "a",tag_suffix = ')') 
  
  return(fin_plot)
}

# Figure 2: Chemistry --------

gen_figure_2a <- function(df){
  
  ### Testing ###
  # df <- CTC_data
  ### Testing ###
  
  df <- df %>% select(Time_AKST,CTC_O3,CTC_NO,CTC_NO2) %>%
    mutate('CTC_NOx' = CTC_NO+CTC_NO2) %>%
    select(-CTC_NO,-CTC_NO2) %>% na.omit() %>% 
    filter(CTC_O3 >= 0)
  
  fin_fig_2a <- ggplot(df,aes(x=CTC_O3,y=CTC_NOx))+
    geom_point()+
    theme_bw()+
    scale_x_continuous(n.breaks = 10,expand = c(0,0))+
    scale_y_continuous(n.breaks = 10,expand = c(0,0))+
    xlab(expression(bold('['*O[3]*'] (ppbv)')))+
    ylab(expression(bold('['*NO[x]*'] (ppbv)')))+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          text = element_text(face = 'bold',family = 'serif',color = 'black'),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10))
  
  return(fin_fig_2a)
}

gen_figure_2b <- function(df,k_val = 100){
  
  ### Testing ###
  # df <- CTC_data
  ### Testing ###
  
  df <- df %>% select(Time_AKST,CTC_O3,CTC_WD,CTC_WS) %>% filter(CTC_O3 > 0)
  
  # Convert the time column to a useable format.
  colnames(df)[colnames(df) == 'Time_AKST'] = 'date'
  colnames(df)[colnames(df) == 'CTC_WS'] = 'ws'
  colnames(df)[colnames(df) == 'CTC_WD'] = 'wd'
  
  font.settings <- list( # for main text (title in this case)
    font = 2, # bold font
    cex = 2, # font magnification based on established fontsize parameter
    fontfamily = "serif") # font
  
  my.theme <- list(
    fontsize = list(text = 30), # Text size for plot labels
    add.line = list(lwd = 2), # Line width for additional lines (like grid lines)
    add.text = list(font = 2, cex = 1, fontfamily = "serif"), # Cardinal directions and other text
    axis.line = list(lwd = 0), # Axis line width (set to 0 to remove scale lines)
    plot.line = list(lwd = 2), # Internal plot lines width
    par.main.text = font.settings, # Font settings for the main title
    axis.text = list(font = 2, cex = 0.7, fontfamily = "serif") # Font settings for axis text (color scale text)
    
    )
  
  fin_fig_2b <- polarPlot(df, pollutant = 'CTC_O3', main = '', width = "fat", auto.text = TRUE,
                 upper = 15, units = 'm s-1', angle.scale = -45, 
                 key.footer = NULL, key.header = expression(bold(O[3]*" (ppbv)")), 
                 exclude.missing = TRUE, 
                 par.settings = my.theme,key.position = 'top',
                 k = k_val)
  dev.off()
  return(fin_fig_2b)
}

gen_figure_2c <- function(df){
  
  ### Testing ###
  # df <- NO3_prod_data
  ### Testing ###
  
  # Load Data ----------
  diel_df <- df %>% select(Time_AKST,rate) %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    mutate_at(vars(Time_AKST),~hour(.))
  
  # Generate Statistics -----------
  
  diel_stat <- diel_df %>% reframe(
    across(rate,calc_stats,.unpack = T)
    ,.by=Time_AKST
  ) %>% pivot_wider(id_cols = Time_AKST,names_from = rate_statistic,values_from = rate_val)
  
  # Plot Diel Profiles -----------
  
  ## Mean +/- 95% CI
  
  fin_fig_2c <- ggplot(diel_stat,aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color = '#473C8B',fill = '#473C8B')+
    geom_line(aes(y=mean),linetype = 'dashed',linewidth = 1,color = '#473C8B')+
    geom_point(aes(y=mean),color = 'black',size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    # scale_y_continuous(labels = label_scientific())+
    xlab('Hour of Day (Local Time)')+
    ylab(expression(bold('P('*NO[3]*') (ppbv '*hr^-1*')')))+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          text = element_text(face = 'bold',family = 'serif',color = 'black'),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10))
  
  return(fin_fig_2c)
}

gen_figure_2d <- function(df){
  
  ### Testing ###
  # df <- NO3_rad_data
  ### Testing ###
  
  # Load Data ----------
  diel_df <- df %>% select(Time_AKST,NO3) %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    mutate_at(vars(Time_AKST),~hour(.))
  
  # Generate Statistics -----------
  
  diel_stat <- diel_df %>% reframe(
    across(NO3,calc_stats,.unpack = T)
    ,.by=Time_AKST
  ) %>% pivot_wider(id_cols = Time_AKST,names_from = NO3_statistic,values_from = NO3_val)
  
  # Plot Diel Profiles -----------
  
  ## Mean +/- 95% CI
  
  fin_fig_2d <- ggplot(diel_stat,aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color = '#473C8B',fill = '#473C8B')+
    geom_line(aes(y=mean),linetype = 'dashed',linewidth = 1,color = '#473C8B')+
    geom_point(aes(y=mean),color = 'black',size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(labels = label_scientific())+
    xlab('Hour of Day (Local Time)')+
    ylab(expression(bold('['*NO[3]*'] (molecules '*cm^-3*')')))+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          text = element_text(face = 'bold',family = 'serif',color = 'black'),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10))
  
  return(fin_fig_2d)
}

# Figure 3: Diel Profiles --------
gen_figure_3 <- function(df){
  
  ### Testing ###
  # df <- hr_data
  ### Testing ###
  
  panel_1 <- c('Furan','Furfural','Methylfuran','Methylfurfural')
  panel_2 <- c('CO','NOx')
  panel_3 <- c('Benzene','C8 Aromatics','Toluene')
  panel_4 <- c('BC','TNRM')
  
  colnames(df) <- gsub('_',' ',colnames(df))
  colnames(df) <- gsub('Time AKST','Time_AKST',colnames(df))
  
  p1_col <- c("#BF3EFF", "#458B00", "#00CDCD", "#FF4040")
  p2_col <- c("#FF7F24","#8B1C62")
  p3_col <- c("#CD1076", "#1874CD", "#CDAD00")
  p4_col <- c("#000000", "#CD6090")
  
  # Load Data ----------
  diel_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    select(Time_AKST, all_of(c(panel_1,panel_2,panel_3,panel_4))) %>%
    mutate_at(vars(Time_AKST),~hour(.)) %>%
    pivot_longer(!Time_AKST,names_to = 'Species',values_to = 'Concentration') %>% 
    na.omit() %>% mutate_at(vars(Concentration),~ifelse(Species == 'CO',.*100,.)) %>% 
    mutate_at(vars(Concentration),~ifelse(Species == 'BC',.*5,.))
  
  # Generate Statistics -----------
  
  diel_stat <- diel_df %>% reframe(
    across(Concentration,calc_stats,.unpack = T)
    ,.by=c(Time_AKST,Species)
  ) %>% pivot_wider(id_cols = c(Time_AKST,Species),names_from = Concentration_statistic,values_from = Concentration_val) %>%
    mutate('panel' = ifelse(Species %in% panel_1,'p1',ifelse(Species %in% panel_2,'p2',ifelse(Species %in% panel_3,'p3','p4'))))
  
  # Plot Diurnal/Weekly Profiles -----------
  
  ## Mean +/- 95% CI
  
  pan_1 <- ggplot(diel_stat %>% filter(panel == 'p1'),aes(x=Time_AKST,color = Species,fill = Species))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab('Furanoids (ppbv)')+
    xlab('Hour of Day (Local Time)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p1_col)+
    scale_fill_manual(values=p1_col)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.25,0.82),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.text = element_text(size = 12))
  
  pan_2 <- ggplot(diel_stat %>% filter(panel == 'p2'),
                  aes(x=Time_AKST,color = factor(Species,levels = c("NOx","CO")),
                      fill = factor(Species,levels = c("NOx","CO"))))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(breaks = seq(0,120,20),minor_breaks = seq(0,120,5),limits = c(20,120),expand = c(0,0),
                       sec.axis = sec_axis(name = 'CO (ppmv)',transform = ~./100,breaks = seq(0,1.2,0.2)))+
    ylab(expression(bold(NO[x]*' (ppbv)')))+
    xlab('Hour of Day (Local Time)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p2_col)+
    scale_fill_manual(values=p2_col)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'inside',
          legend.position.inside = c(0.25,0.8),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.text = element_text(size = 12))
  
  pan_3 <- ggplot(diel_stat %>% filter(panel == 'p3'),aes(x=Time_AKST,color = Species,fill = Species))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(breaks = seq(0,16,2),minor_breaks = seq(0,16,0.5),limits = c(0,15),expand = c(0,0))+
    xlab('Hour of Day (Local Time)')+
    ylab('Aromatics (ppbv)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p3_col)+
    scale_fill_manual(values=p3_col)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          legend.position = 'inside',legend.position.inside = c(0.25,0.8),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.text = element_text(size = 12))
  
  pan_4 <- ggplot(diel_stat %>% filter(panel == 'p4'),aes(x=Time_AKST,color = Species,fill = Species))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,30,5),minor_breaks = seq(0,30,1),limits = c(0,26),
                       sec.axis = sec_axis(name = expression(bold('BC ('*mu*'g '*m^-3*')')),
                                           transform = ~./5,breaks = seq(0,6,1)))+
    xlab('Hour of Day (Local Time)')+
    ylab(expression(bold('TNRM ('*mu*'g '*m^-3*')')))+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p4_col)+
    scale_fill_manual(values=p4_col)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          legend.position = 'inside',legend.position.inside = c(0.25,0.8),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.text = element_text(size = 12))
  
  fin_fig_3 <- ((pan_1+pan_2)/(pan_3+pan_4) & xlab(NULL) & theme(plot.tag.position  = c(0.15, .96)))+
    plot_annotation(tag_levels = "a",tag_suffix = ')') 
  
  fin_fig_3 <- wrap_elements(panel = fin_fig_3) +
    labs(tag = "Hour of Day (Local Time)") +
    theme(plot.tag = element_text(size = rel(1),face = 'bold',family = 'serif'),
          plot.tag.position = "bottom")
  return(fin_fig_3)
}

# Figure 4: Scatter Plots ---------

gen_figure_4 <- function(df){
  
  ### Testing ###
  # df <- hr_data
  ### Testing ###
  
  df <- df %>% mutate('sum_aromatics' = Benzene+Toluene+C8_Aromatics+C9_Aromatics+Naphthalene+C10_Aromatics)
  
  panel_1 <- c('Benzene','Toluene')
  panel_2 <- c('Furan','Methylfuran')
  panel_3 <- c('sum_aromatics','TVOC')
  panel_4 <- c('TVOC','CO')
  
  geom_eq <- function(x,y,coord_x,coord_y){
    
    reg_model <- lm(y~x)
    reg_model_summary <- summary(reg_model)
    reg_stats <- c(as.numeric(reg_model_summary$coefficients),as.numeric(reg_model_summary$adj.r.squared))
    if(reg_stats[1]<0){
      reg_eq <- paste0(
        "y == ", signif(reg_stats[2],3), " * x ", round(reg_stats[1], 2), 
        " ~ italic(R)^2 == ", round(reg_stats[length(reg_stats)], 2)
      )
    }else{
      reg_eq <- paste0(
        "y == ", signif(reg_stats[2],3), " * x + ", round(reg_stats[1], 2), 
        " ~ italic(R)^2 == ", round(reg_stats[length(reg_stats)], 2)
      )
    }
    reg_line <- geom_abline(intercept = reg_stats[1],slope = reg_stats[2],color='blue',linewidth = 2)
    reg_label <- annotate('text',x=coord_x,y=coord_y,label = reg_eq,parse = T,size = 5)
    return(list(reg_line,reg_label))
  }
  
  p1 <- ggplot(df %>% select(all_of(panel_1)),aes(x=Benzene,y=Toluene))+
    geom_point(size = 2)+
    geom_eq(x=df$Benzene,y=df$Toluene,coord_x = 3,coord_y = 30)+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'))+
    xlab('Benzene (ppbv)')+
    ylab('Toluene (ppbv)')
  
  p2 <- ggplot(df %>% select(all_of(panel_2)),aes(x=Furan,y=Methylfuran))+
    geom_point(size = 2)+
    geom_eq(x=df$Furan,y=df$Methylfuran,coord_x = 0.12,coord_y = 0.4)+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'))+
    xlab('Furan (ppbv)')+
    ylab('Methylfuran (ppbv)')
  
  p3 <- ggplot(df %>% select(all_of(panel_3)),aes(x=TVOC,y=sum_aromatics))+
    geom_point(size = 2)+
    geom_eq(x=df$TVOC,y=df$sum_aromatics,coord_x = 60,coord_y = 105)+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'))+
    xlab('TVOC (ppbv)')+
    ylab(expression(bold(Sigma*'Aromatics (ppbv)')))
  
  p4 <- ggplot(df %>% select(all_of(panel_4)),aes(x=CO,y=TVOC))+
    geom_point(size = 2)+
    geom_eq(x=df$CO,y=df$TVOC,coord_x = 2,coord_y = 20)+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'))+
    xlab('CO (ppmv)')+
    ylab('TVOC (ppbv)')
  
 fin_fig_4 <- ((p3+p4)/(p1+p2) & theme(plot.tag.position  = c(0.15, .96)))+
   plot_annotation(tag_levels = "a",tag_suffix = ')')
 
  return(fin_fig_4)
}

# Figure 5: Pie Chart ---------
gen_figure_5 <- function(df,met_df){
  
  ### Testing ###
  # df <- hr_data
  # met_df <- met_data
  ### Testing ###
  
  ## Select Species 
  spec_list <- c('Time_AKST','Methanol','Acetonitrile','Acetaldehyde','Formic_Acid',
                 'Ethanol','Methanethiol','Acetone','Acetic_Acid','Furan','Isoprene',
                 'MVK_MACR','Methyl_Glyoxal','MEK','Hydroxyacetone','Benzene','Methylfuran',
                 'Pentanone','Toluene','Furfural','Maleic_Anhydride','Hexanone','C8_Aromatics'
                 ,'Methylfurfural','C9_Aromatics','Naphthalene','C10_Aromatics','Monoterpenes'
                 ,'D5_Siloxane','Formaldehyde','CO','CO2','CH4','Org')
  
  MW_data <- data.frame('Species' = c('Methanol','Acetonitrile','Acetaldehyde','Formic_Acid',
                        'Ethanol','Methanethiol','Acetone','Acetic_Acid','Furan','Isoprene',
                        'MVK_MACR','Methyl_Glyoxal','MEK','Hydroxyacetone','Benzene','Methylfuran',
                        'Pentanone','Toluene','Furfural','Maleic_Anhydride','Hexanone','C8_Aromatics',
                        'Methylfurfural','C9_Aromatics','Naphthalene','C10_Aromatics','Monoterpenes',
                        'D5_Siloxane','Formaldehyde','CO','CO2','CH4','Org')
                        ,'MW' = c(32.042,41.0527,44.053,46.025,46.069,48.108,48.108,60.052,68.075,
                                  68.119,70.091,72.063,72.107,74.079,78.114,82.102,86.134,92.141,
                                  96.085,98.057,100.161,106.168,110.112,120.195,128.174,134.222,
                                  136.238,370.775,30.011,28.01,44.01,16.04,NA)
                        ,'Class' = c('Alcohols','Other NMVOCs','Other NMVOCs','Other NMVOCs',
                                     'Alcohols','Other NMVOCs','Other NMVOCs','Other NMVOCs','Other NMVOCs','Other NMVOCs',
                                     'Other NMVOCs','Other NMVOCs','Other NMVOCs','Other NMVOCs','Aromatics','Other NMVOCs',
                                     'Other NMVOCs','Aromatics','Other NMVOCs','Other NMVOCs','Other NMVOCs','Aromatics',
                                     'Other NMVOCs','Aromatics','Aromatics','Aromatics','Other NMVOCs',
                                     'Other NMVOCs','Other NMVOCs','CO','CO2','CH4','OA'))
  ## Convert from ppbv or ppmv to non-scaled mixing ratio
  df <- df %>% select(all_of(spec_list)) %>% 
    mutate_at(vars(-Time_AKST,-CO,-CO2,-CH4,-Org),~./(10^9)) %>%
    mutate_at(vars(CO,CO2,CH4),~./(10^6))
  
  ## Merge meteorological data and molecular weight data to mixing ratio data
  df <- df %>% pivot_longer(cols = !Time_AKST,names_to = 'Species',values_to = 'Concentration') %>%
    merge(.,MW_data,by = 'Species',all.x=T) %>%
    merge(.,met_data %>% select(Time_AKST,altimeter_set_1,air_temp_set_1),by = 'Time_AKST',all.x=T)
  
  colnames(df)[colnames(df) == 'altimeter_set_1'] <- 'P_Pa'
  colnames(df)[colnames(df) == 'air_temp_set_1'] <- 'T_C'
  
  ## Convert from Mixing Ratio to Mass Density
  df_px <- px_Cx_df(df)
  ## Calculate Statistics
  df_stat <- df_px %>% reframe(
    across(Concentration,calc_stats,.unpack = T)
    ,.by=Species
  ) %>% filter(Concentration_statistic == 'mean') %>% select(-Concentration_statistic)
  colnames(df_stat)[colnames(df_stat) == 'Concentration_val'] <- 'mean'
  VOC_PM_ROC <- sum((df_stat %>% filter(!Species %in% c('CO','CO2','CH4')))$mean)
  
  ## Create Total ROC and Percent ROC by Mass Columns
  df_stat <- df_stat  %>% filter(Species != 'CO2') %>%
    mutate('Phase' = ifelse(Species %in% c('Org','BC'),'Aerosol',ifelse(Species %in% c('CO','CO2','CH4'),Species,'NMVOC'))) %>% 
    mutate('vpROC' = VOC_PM_ROC) %>%
    mutate('vpROC_perc' = ifelse(Phase %in% c('CO','CO2','CH4'),NA,(mean/vpROC)*100)) %>%
    merge(.,MW_data,by='Species',all.x = T)
  
  ## Plot ROC Pie Charts
  ### NMVOCs and Aerosols (OA+BC)
  vpROC_df <- df_stat %>% filter(Phase %in% c('Aerosol','NMVOC')) %>%
    select(Class,Phase,vpROC,vpROC_perc) %>% reframe(
      across(!Phase,sum,.unpack = T)
      ,.by=Class
    )
  
  tot_ROC <- round(df_stat$vpROC[[1]],2)
  
  piecol <- c("#6B8E23", "#4F94CD", "#FFC125", "#7A67EE")
  
  pie_1 <- ggplot(vpROC_df, aes(x = "", y = vpROC_perc, fill = Class))+
    geom_col()+
    geom_text(aes(label = paste0(round(vpROC_perc,0),'%')),
              position = position_stack(vjust = 0.5),
              fontface = 'bold',color = 'black',size = 8) +
    coord_polar(theta = "y")+ 
    theme_void()+
    theme(legend.title = element_blank(),
          legend.text = element_text(face = 'bold',family = 'serif',
                                     color='black',size = 11))+
    scale_fill_manual(values = piecol)
  
  return(pie_1)
}

# Figure 6: PMF Diurnals and Bar Chart ----------

gen_figure_6a <- function(df){
  
  ### Testing ###
  # df <- PMF_TS_data
  ### Testing ###
  
  PMF_factors <- c('BB-SL')
  
  factor_col <- c("#FF7F24")
  
  
  diel_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    select(Time_AKST, all_of(PMF_factors)) %>%
    mutate_at(vars(Time_AKST),~hour(.)) %>%
    pivot_longer(!Time_AKST,names_to = 'Factors',values_to = 'Contributions') %>% 
    na.omit()
  
  diel_stat <- diel_df %>% reframe(
    across(Contributions,calc_stats,.unpack = T)
    ,.by=c(Time_AKST,Factors)
  ) %>% pivot_wider(id_cols = c(Time_AKST,Factors),
                    names_from = Contributions_statistic,
                    values_from = Contributions_val)
  
  
  
  ## Mean +/- 95% CI
  
  fin_fig_6a <- ggplot(diel_stat,aes(x=Time_AKST,color = Factors,fill = Factors))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1)+
    geom_line(aes(y=mean),linetype = 'dashed',linewidth = 1)+
    geom_point(aes(y=mean),color = 'black',size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),
                       minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(minor_breaks = seq(0,3,0.1))+
    xlab('Hour of Day (Local Time)')+
    ylab('Normalized Contributions (mean = 1)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=factor_col)+
    scale_fill_manual(values=factor_col)+
    facet_wrap(~Factors,scales='free_y')+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE),
           color = guide_legend(nrow = 1))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = 'inside',
          legend.position.inside = c(0.2,0.85),
          legend.background = element_blank())
  
  return(fin_fig_6a)
}

gen_figure_6b <- function(df){
  
  ### Testing ###
  # df <- PMF_TS_data
  ### Testing ###
  
  PMF_factors <- c('Traffic')
  
  factor_col <- c("#436EEE")
  
  
  diel_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    select(Time_AKST, all_of(PMF_factors)) %>%
    mutate_at(vars(Time_AKST),~hour(.)) %>%
    pivot_longer(!Time_AKST,names_to = 'Factors',values_to = 'Contributions') %>% 
    na.omit()
  
  diel_stat <- diel_df %>% reframe(
    across(Contributions,calc_stats,.unpack = T)
    ,.by=c(Time_AKST,Factors)
  ) %>% pivot_wider(id_cols = c(Time_AKST,Factors),names_from = Contributions_statistic,values_from = Contributions_val)
  
  
  
  ## Mean +/- 95% CI
  
  fin_fig_6b <- ggplot(diel_stat,aes(x=Time_AKST,color = Factors,fill = Factors))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1)+
    geom_line(aes(y=mean),linetype = 'dashed',linewidth = 1)+
    geom_point(aes(y=mean),color = 'black',size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(minor_breaks = seq(0,3,0.1))+
    xlab('Hour of Day (Local Time)')+
    ylab('Normalized Contributions (mean = 1)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=factor_col)+
    scale_fill_manual(values=factor_col)+
    facet_wrap(~Factors,scales='free_y')+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE),
           color = guide_legend(nrow = 1))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = 'inside',
          legend.position.inside = c(0.2,0.85),
          legend.background = element_blank())
  
  return(fin_fig_6b)
}

gen_figure_6c <- function(df){
  
  ### Testing ###
  # df <- PMF_TS_data
  ### Testing ###
  
  PMF_factors <- c('BB-LL')
  
  factor_col <- c("#EE0000")
  
  
  ts_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    select(Time_AKST, all_of(PMF_factors)) %>% ts_gap() %>%
    pivot_longer(!Time_AKST,names_to = 'Factors',values_to = 'Contributions')
  
  ## Mean +/- 95% CI
  
  fin_fig_6c <- ggplot(ts_df,aes(x=Time_AKST,color = Factors))+
    geom_line(aes(y=Contributions),linewidth = 1)+
    theme_bw()+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day')+
    scale_y_continuous(limits = c(0,NA),expand=c(0,0))+
    xlab('Time (AKST)')+
    ylab('Normalized g-values (mean = 1)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=factor_col)+
    scale_fill_manual(values=factor_col)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE),
           color = guide_legend(nrow = 1))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          legend.title = element_blank(),
          legend.position = 'inside',
          legend.position.inside = c(0.8625,0.8),
          legend.background = element_blank(),
          legend.key=element_blank(),
          axis.title.y = element_blank())
  
  return(fin_fig_6c)
}

gen_figure_6d <- function(df){
  
  ### Testing ###
  # df <- PMF_TS_data
  ### Testing ###
  
  PMF_factors <- c('HO')
  
  factor_col <- c("#8B2500")
  
  
  ts_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    select(Time_AKST, all_of(PMF_factors)) %>% ts_gap() %>%
    pivot_longer(!Time_AKST,names_to = 'Factors',values_to = 'Contributions')
  
  ## Mean +/- 95% CI
  
  fin_fig_6d <- ggplot(ts_df,aes(x=Time_AKST,color = Factors))+
    geom_line(aes(y=Contributions),linewidth = 1)+
    theme_bw()+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day')+
    scale_y_continuous(limits = c(0,NA),expand=c(0,0))+
    xlab('Time (AKST)')+
    ylab('Normalized g-values (mean = 1)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=factor_col)+
    scale_fill_manual(values=factor_col)+
    guides(y = guide_axis(minor.ticks = TRUE),
           x = guide_axis(minor.ticks = TRUE),
           color = guide_legend(nrow = 1))+
    theme(axis.ticks.length = unit(5, "pt"),axis.minor.ticks.length = rel(0.5),
          legend.title = element_blank(),
          legend.position = 'inside',
          legend.position.inside = c(0.85,0.8),
          legend.background = element_blank(),
          legend.key=element_blank(),
          axis.title.y = element_blank())
  
  return(fin_fig_6d)
}

# Figure 7: PMF Factor Contributions to Species Total Fit -----------

gen_figure_7 <- function(df){
  
  ### Testing ###
  # df <- tabs4
  ### Testing ###
  
  PMF_factors <- c('Aged','BB-LL','BB-SL','HO',expression(bold(O[3])),'Traffic')
  
  factor_col <- c('#556B2F',"#EE0000","#FF7F24","#8B2500",'#8B1C62',"#436EEE")
  
  spec_list <- c('TNRM','SO4','NOx','CO','Formaldehyde','Formic_Acid',
                 'Acetic_Acid','Furfural','Furan','Benzene','Toluene')
  
  df <- df %>% filter(Species_Name %in% spec_list)
  
  df <- df %>% mutate_at(vars(Species_Name),~gsub('_',' ',.))
  
  spec_labs <- c(expression(bold('Acetic Acid')),expression(bold('Benzene')),
                 expression(bold('CO')),expression(bold('Formaldehyde')),
                 expression(bold('Formic Acid')),expression(bold('Furan')),
                 expression(bold('Furfural')),expression(bold(NO[x])),
                 expression(bold({SO[4]^{'2-'}~Aerosol})),expression(bold('TNRM - Aerosol')),
                 expression(bold('Toluene')))
  
  fin_fig_7 <- ggplot(df,aes(x=Species_Name,y=perc_cont))+
    geom_col(aes(fill = PMF_factor),color = 'black')+
    scale_y_continuous(breaks = seq(0,100,10),
                       ,minor_breaks = seq(0,100,5),expand = c(0,0))+
    scale_fill_manual(values = factor_col,labels = PMF_factors)+
    scale_x_discrete(labels = spec_labs)+
    xlab('')+
    ylab('% of Species Sum')+
    theme_bw()+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'),
          legend.position = 'top',legend.title = element_blank(),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5,size = 12),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 15))+
    guides(fill = guide_legend(nrow=1))
  
  return(fin_fig_7)
}

# Supplemental Information ---------------------

# Figure S1: Anticorrelations -----------

gen_figure_1_SI <- function(df){
  
  ### Testing ###
  # df <- CTC_data
  ### Testing ###
  
  df <- df %>% select(Time_AKST,CTC_O3,CTC_NO,CTC_NO2,PM1,CTC_CO) %>%
    mutate('CTC_NOx' = CTC_NO+CTC_NO2) %>%
    select(-CTC_NO,-CTC_NO2) %>% na.omit() %>% 
    filter(CTC_O3 >= 0)
  
  pan_a <- ggplot(df,aes(x=CTC_O3,y=CTC_CO))+
    geom_point()+
    theme_bw()+
    scale_x_continuous(name = expression(bold('['*O[3]*'] (ppbv)')),
                       n.breaks = 10,expand = c(0,0))+
    scale_y_continuous(name = '[CO] (ppmv)',
                       n.breaks = 10,expand = c(0,0),limits = c(0,5))+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10))
  
  pan_b <- ggplot(df,aes(x=CTC_O3,y=PM1))+
    geom_point()+
    theme_bw()+
    scale_x_continuous(name = expression(bold('['*O[3]*'] (ppbv)')),
                       n.breaks = 10,expand = c(0,0))+
    scale_y_continuous(name = expression(bold(PM[1]*' ('*mu*'g '*m^-3*')')),
                       n.breaks = 10,expand = c(0,0))+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10))
  
  fin_ac_SI <- ((pan_a + pan_b) & theme(plot.tag.position  = c(0.15, .96)))+
    plot_annotation(tag_levels = "a",tag_suffix = ')') 
  
  return(fin_ac_SI)
}

# Figure S2: Diel Profile of Toluene/Benzene Ratio ----------

gen_figure_2_SI <- function(df){
  
  ### Testing ###
  # df <- hr_data
  ### Testing ###
  
  # Load Data ----------
  diel_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    mutate('TBR' = Toluene/Benzene) %>%
    select(Time_AKST, TBR) %>%
    mutate_at(vars(Time_AKST),~hour(.)) %>%
    pivot_longer(!Time_AKST,names_to = 'Species',values_to = 'Concentration') %>% 
    na.omit()
  
  # Generate Statistics -----------
  
  diel_stat <- diel_df %>% reframe(
    across(Concentration,calc_stats,.unpack = T)
    ,.by=c(Time_AKST,Species)
  ) %>% pivot_wider(id_cols = c(Time_AKST,Species),
                    names_from = Concentration_statistic,
                    values_from = Concentration_val)
  
  # Plot Diurnal/Weekly Profiles -----------
  
  ## Mean +/- 95% CI
  
  fin_TBR_fig <- ggplot(diel_stat,aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),
                alpha = 0.5,linewidth = 1,color = '#cc5e3a',fill = '#cc5e3a')+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,
              color = 'black')+
    geom_point(aes(y=mean),size = 2,color = 'black')+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    scale_y_continuous(breaks = seq(1.5,2.6,0.1),minor_breaks = seq(1.5,2.6,0.02),limits = c(1.5,2.6),expand = c(0,0))+
    ylab(expression(bold('Toluene/Benzene (ppbv '*ppbv^-1*')')))+
    xlab('Hour of Day (Local Time)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10))
  
  return(fin_TBR_fig)
}

# Figure S3: ADEC vs PMF Factor Source Profiles for RWC and HO -----------

gen_figure_3_SI <- function(df){
  
  ### Testing ###
  # df <- SI_plot_data
  ### Testing ###
  
  panel_a <- c('ADEC_HO_103','PMF_HO_6F')
  panel_b <- c('ADEC_RWC_107','PMF_BB_6F')
  panel_c <- c('ADEC_HO_103','PMF_HO_7F')
  panel_d <- c('ADEC_RWC_107','PMF_BB_7F')
  panel_e <- c('ADEC_HO_103','PMF_HO1_9F')
  panel_f <- c('ADEC_HO_103','PMF_HO2_9F')
  
  geom_eq <- function(x,y,coord_x,coord_y){
    
    reg_model <- lm(y~x)
    reg_model_summary <- summary(reg_model)
    reg_stats <- c(as.numeric(reg_model_summary$coefficients),as.numeric(reg_model_summary$adj.r.squared))
    if(reg_stats[1]<0){
      reg_eq <- bquote(y == .(signif(reg_stats[2], 3)) * x * .(round(reg_stats[1], 2)) * ";" ~ italic(R)^2 == .(round(reg_stats[length(reg_stats)], 2)))
    } else{
      reg_eq <- bquote(y == .(signif(reg_stats[2], 3)) * x + .(round(reg_stats[1], 2)) * ";" ~ italic(R)^2 == .(round(reg_stats[length(reg_stats)], 2)))
    }
    reg_line <- geom_abline(intercept = reg_stats[1],slope = reg_stats[2],color='blue',linewidth = 2)
    reg_label <- annotate('text',x=coord_x,y=coord_y,label = reg_eq,parse = F,size = 4)
    return(list(reg_line,reg_label))
  }
  
  p1 <- ggplot(df %>% select(all_of(panel_a)),aes(x=ADEC_HO_103,y=PMF_HO_6F))+
    geom_point(size = 2)+
    geom_eq(x=df$ADEC_HO_103,y=df$PMF_HO_6F,coord_x = 0.4,coord_y = 0.2)+
    theme_bw()+
    annotate('text',x=0.025,y=1.85,label = 'a)',size = 5,fontface = 'bold')+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    ylab('6F PMF: HO')+
    theme(axis.title.y.left = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'),
          axis.title.y.right = element_blank(),
          axis.title.x = element_blank())
  
  p2 <- ggplot(df %>% select(all_of(panel_b)),aes(x=ADEC_RWC_107,y=PMF_BB_6F))+
    geom_point(size = 2)+
    geom_eq(x=df$ADEC_RWC_107,y=df$PMF_BB_6F,coord_x = 0.6,coord_y = 0.075)+
    theme_bw()+
    annotate('text',x=0.025,y=0.665,label = 'b)',size = 5,fontface = 'bold')+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10,position = "right")+
    theme(axis.title.y.right = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'),
          axis.title.y.left = element_blank(),
          axis.title.x = element_blank())+
    ylab('6F PMF: BB-SL + BB-LL')
  
  p3 <- ggplot(df %>% select(all_of(panel_c)),aes(x=ADEC_HO_103,y=PMF_HO_7F))+
    geom_point(size = 2)+
    geom_eq(x=df$ADEC_HO_103,y=df$PMF_HO_7F,coord_x = 0.4,coord_y = 0.2)+
    theme_bw()+
    annotate('text',x=0.025,y=1.6,label = 'c)',size = 5,fontface = 'bold')+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    theme(axis.title.y.left = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'),
          axis.title.y.right = element_blank(),
          axis.title.x = element_blank())+
    ylab('7F PMF: HO')
  
  p4 <- ggplot(df %>% select(all_of(panel_d)),aes(x=ADEC_RWC_107,y=PMF_BB_7F))+
    geom_point(size = 2)+
    geom_eq(x=df$ADEC_RWC_107,y=df$PMF_BB_7F,coord_x = 0.6,coord_y = 0.075)+
    theme_bw()+
    annotate('text',x=0.025,y=0.585,label = 'd)',size = 5,fontface = 'bold')+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10,position = "right")+
    theme(axis.title.y.right = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'),
          axis.title.y.left = element_blank(),
          axis.title.x = element_blank())+
    ylab('7F PMF: BB-SL + BB-LL')
  
  p5 <- ggplot(df %>% select(all_of(panel_e)),aes(x=ADEC_HO_103,y=PMF_HO1_9F))+
    geom_point(size = 2)+
    geom_eq(x=df$ADEC_HO_103,y=df$PMF_HO1_9F,coord_x = 0.4,coord_y = 0.075)+
    theme_bw()+
    annotate('text',x=0.025,y=0.8,label = 'e)',size = 5,fontface = 'bold')+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10)+
    theme(axis.title.y.left = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.title.x = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.title.y.right = element_blank(),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'))+
    xlab('ADEC HO #103')+
    ylab('9F PMF: HO-1')
  
  p6 <- ggplot(df %>% select(all_of(panel_f)),aes(x=ADEC_HO_103,y=PMF_HO2_9F))+
    geom_point(size = 2)+
    geom_eq(x=df$ADEC_HO_103,y=df$PMF_HO2_9F,coord_x = 0.45,coord_y = 0.2)+
    theme_bw()+
    annotate('text',x=0.025,y=1.85,label = 'f)',size = 5,fontface = 'bold')+
    scale_x_continuous(n.breaks = 10)+
    scale_y_continuous(n.breaks = 10,position = "right")+
    theme(axis.title.y.right = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.title.x = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.title.y.left = element_blank(),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          text = element_text(family = 'serif',face = 'bold',colour = 'black'))+
    xlab('ADEC RWC #107 (b,d); ADEC HO #103 (f)')+
    ylab('9F PMF: HO-2')
  
  fin_fig_S3 <- ((p1+p2)/(p3+p4)/(p5+p6)) #& theme(plot.tag.position  = c(0.15, .96)))+
    #plot_annotation(tag_levels = "a",tag_suffix = ')')
  
  return(fin_fig_S3)
}
  
# Figure S4: Sulfate Diurnal Profile by Month -------------

gen_figure_4_SI <- function(df){
  
  ### Testing ###
  # df <- hr_data
  ### Testing ###
  
  # Load Data ----------
  diel_df <- df %>% mutate('mon' = lubridate::month(Time_AKST,label = T,abbr = F)) %>%
    mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    select(Time_AKST,mon,SO4) %>%
    mutate_at(vars(Time_AKST),~lubridate::hour(.)) %>%
    pivot_longer(!c(Time_AKST,mon),names_to = 'Species',values_to = 'Concentration') %>% 
    na.omit()
  
  # Generate Statistics -----------
  
  diel_stat <- diel_df %>% reframe(
    across(Concentration,calc_stats,.unpack = T)
    ,.by=c(Time_AKST,mon,Species)
  ) %>% pivot_wider(id_cols = c(Time_AKST,mon,Species),
                    names_from = Concentration_statistic,
                    values_from = Concentration_val) %>%
    mutate('plot_lab' = ifelse(mon == 'January','a)','b)'))
  
  # Plot Diurnal/Weekly Profiles -----------
  
  ## Mean +/- 95% CI
  
  fin_SO4_fig <- ggplot(data = diel_stat,aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),
                alpha = 0.5,linewidth = 1,color = '#EEB422',fill = '#EEB422')+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,
              color = 'black')+
    geom_point(aes(y=mean),size = 2,color = 'black')+
    theme_bw()+
    geom_text(aes(x=2,y=6.5,label = plot_lab),color = "black",fontface = "bold",size = 4)+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0.003,0))+
    scale_y_continuous(breaks = seq(0,7,1),minor_breaks = seq(0,7,0.1),limits = c(0,7),expand = c(0.003,0.003))+
    ylab(expression(bold({SO[4]^{'2-'}})))+
    xlab('Hour of Day (Local Time)')+
    facet_wrap(~mon)+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold',family = 'serif',color = 'black',size = 12))
  
  return(fin_SO4_fig)
}

# Figure S5: PMF Time Series and Diel Profiles --------
gen_figure_5_SI <- function(df){
  
  ### Testing ###
  # df <- PMF_TS_data
  ### Testing ###
  
  panel_1 <- 'Aged'
  panel_2 <- 'BB-LL'
  panel_3 <- 'BB-SL'
  panel_4 <- 'HO'
  panel_5 <- 'O3'
  panel_6 <- 'Traffic'
  
  p1_col <- '#556B2F'
  p2_col <- "#EE0000"
  p3_col <- "#FF7F24"
  p4_col <- "#8B2500"
  p5_col <- '#8B1C62'
  p6_col <- "#436EEE"
  
  # Load Data ----------
  diel_df <- df %>% mutate_at(vars(Time_AKST),~floor_date(.,unit = 'hour')) %>%
    mutate_at(vars(Time_AKST),~hour(.)) %>%
    pivot_longer(!Time_AKST,names_to = 'Factor',values_to = 'g_value')
  
  # Generate Statistics -----------
  
  diel_stat <- diel_df %>% reframe(
    across(g_value,calc_stats,.unpack = T)
    ,.by=c(Time_AKST,Factor)
  ) %>% pivot_wider(id_cols = c(Time_AKST,Factor),names_from = g_value_statistic,values_from = g_value_val)
  
  # Plot Diurnal/Weekly Profiles -----------
  
  ## Mean +/- 95% CI
  
  pan_1 <- ggplot(diel_stat %>% filter(Factor == panel_1),aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color= p1_col,fill= p1_col)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,color= p1_col)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab(NULL)+
    xlab(NULL)+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p1_col)+
    scale_fill_manual(values=p1_col)+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.25,0.82),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.text = element_text(size = 12))
  
  pan_2 <- ggplot(diel_stat %>% filter(Factor == panel_2),aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color= p2_col,fill= p2_col)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,color= p2_col)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab(NULL)+
    xlab(NULL)+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p2_col)+
    scale_fill_manual(values=p2_col)+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'none')
  
  pan_3 <- ggplot(diel_stat %>% filter(Factor == panel_3),aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color= p3_col,fill= p3_col)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,color= p3_col)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab(NULL)+
    xlab(NULL)+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p3_col)+
    scale_fill_manual(values=p3_col)+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'none')
  
  pan_4 <- ggplot(diel_stat %>% filter(Factor == panel_4),aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color= p4_col,fill= p4_col)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,color= p4_col)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab(NULL)+
    xlab(NULL)+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p4_col)+
    scale_fill_manual(values=p4_col)+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'none')
  
  pan_5 <- ggplot(diel_stat %>% filter(Factor == panel_5),aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color= p5_col,fill= p5_col)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,color= p5_col)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab(NULL)+
    xlab(NULL)+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p5_col)+
    scale_fill_manual(values=p5_col)+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          legend.position = 'none')
  
  pan_6 <- ggplot(diel_stat %>% filter(Factor == panel_6),aes(x=Time_AKST))+
    geom_ribbon(aes(ymin = mean-CI95_ME,ymax = mean+CI95_ME),alpha = 0.5,linewidth = 1,color= p6_col,fill= p6_col)+
    geom_line(aes(y=mean),linetype = 'solid',linewidth = 1,color= p6_col)+
    geom_point(aes(y=mean),size = 2)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0,24,4),labels = seq(0,24,4),minor_breaks = seq(0,24,1),expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0,0.2,0.05),minor_breaks = seq(0,0.2,0.01),limits = c(0,0.2),expand = c(0,0))+
    ylab(NULL)+
    xlab('Hour of Day (Local Time)')+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'))+
    scale_color_manual(values=p6_col)+
    scale_fill_manual(values=p6_col)+
    theme(axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title.x = element_text(size = rel(1),face = 'bold',family = 'serif',angle=0),
          legend.position = 'none')
  
  fin_diel <- ((pan_1 & theme(plot.tag.position  = c(0.16, .94)))/
                 (pan_2 & theme(plot.tag.position  = c(0.16, .94)))/
                 (pan_3 & theme(plot.tag.position  = c(0.16, .94)))/
                 (pan_4 & theme(plot.tag.position  = c(0.16, .94)))/
                 (pan_5 & theme(plot.tag.position  = c(0.16, .94)))/
                 (pan_6 & theme(plot.tag.position  = c(0.16, .94))) & ylab(NULL))+
    plot_annotation(tag_levels = "a",tag_suffix = ')')
  
  fin_diel <- wrap_elements(panel = fin_diel) +
    labs(tag = expression(bold("Normalized g-values ("*bar(x)*' = 1)'))) +
    #xlab('Hour of Day (Local Time)')+
    theme(plot.tag = element_text(size = rel(1),face = 'bold',family = 'serif',angle=90),
          plot.tag.position = "left")+plot_layout(axes = "collect")
  
  time_seq <- data.frame('Time_AKST' = seq.POSIXt(from = as.POSIXct('2022-01-19 00:00:00',tz='America/Sitka'),
                         as.POSIXct('2022-02-25 12:00:00',tz='America/Sitka'),by = '1 hour'))
  
  df <- merge(df,time_seq,by = 'Time_AKST',all = T)
  
  ts_1 <- ggplot(df,aes(x=Time_AKST,y=Aged,color = p1_col))+
    geom_line(linewidth=1)+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    ylab(NULL)+
    scale_y_continuous(limits = c(0,4),expand = c(0,0))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.90,0.90),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p1_col,labels = panel_1)+
    guides(color = guide_legend(nrow = 1))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  ts_2 <- ggplot(df,aes(x=Time_AKST,y=`BB-LL`,color = p2_col))+
    geom_line(linewidth=1)+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    ylab(NULL)+
    scale_y_continuous(limits = c(0,3.5),expand = c(0,0))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.90,0.90),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p2_col,labels = panel_2)+
    guides(color = guide_legend(nrow = 1))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  ts_3 <- ggplot(df,aes(x=Time_AKST,y=`BB-SL`,color = p3_col))+
    geom_line(linewidth=1)+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    ylab(NULL)+
    scale_y_continuous(limits = c(0,7.5),expand = c(0,0))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.90,0.90),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p3_col,labels = panel_3)+
    guides(color = guide_legend(nrow = 1))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  ts_4 <- ggplot(df,aes(x=Time_AKST,y=HO,color = p4_col))+
    geom_line(linewidth=1)+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    ylab(NULL)+
    scale_y_continuous(limits = c(0,5),expand = c(0,0))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.90,0.90),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p4_col,labels = panel_4)+
    guides(color = guide_legend(nrow = 1))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  ts_5 <- ggplot(df,aes(x=Time_AKST,y=O3,color = p5_col))+
    geom_line(linewidth=1)+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    ylab(NULL)+
    scale_y_continuous(limits = c(0,3.5),expand = c(0,0))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          legend.position = 'inside',legend.position.inside = c(0.90,0.90),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p5_col,labels = panel_5)+
    guides(color = guide_legend(nrow = 1))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  ts_6 <- ggplot(df,aes(x=Time_AKST,y=Traffic,color = p6_col))+
    geom_line(linewidth=1)+
    scale_x_datetime(date_breaks = '1 week',date_minor_breaks = '1 day',expand = c(0,0))+
    theme_bw()+
    ylab(NULL)+
    xlab('Time (AKST)')+
    scale_y_continuous(limits = c(0,13),breaks = seq(0,14,by=2),minor_breaks = seq(0,14,by=0.5),expand = c(0,0))+
    theme(legend.position = 'inside',legend.position.inside = c(0.90,0.90),
          legend.title = element_blank(),legend.background = element_blank(),
          legend.key=element_blank())+
    scale_color_manual(values = p6_col,labels = panel_6)+
    guides(color = guide_legend(nrow = 1))+
    theme(text = element_text(family = 'serif',face = 'bold',color = 'black'),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12))
  
  fin_ts <- ((ts_1 & theme(plot.tag.position  = c(0.077, .93)))/
               (ts_2 & theme(plot.tag.position  = c(0.077, .93)))/
               (ts_3 & theme(plot.tag.position  = c(0.077, .93)))/
               (ts_4 & theme(plot.tag.position  = c(0.077, .93)))/
               (ts_5 & theme(plot.tag.position  = c(0.077, .93)))/
               (ts_6 & theme(plot.tag.position  = c(0.077, .93))) & ylab(NULL))+
    plot_annotation(tag_levels = list(c('g','h','i','j','k','l'),1),tag_suffix = ')')
  
  fin_ts <- wrap_elements(panel = fin_ts)
  
  fin_fig_S5 <- fin_diel+fin_ts+plot_layout(widths = c(3,10))
  
  return(fin_fig_S5)
}

# Figure S6: PMF Factor Profiles (Fingerprints) --------
gen_figure_6_SI <- function(df){
  
  ### Testing ###
  # df <- tabs4
  ### Testing ###
  
  PMF_factors <- c('Aged','BB-LL','BB-SL','HO',expression(bold(O[3])),'Traffic')
  
  factor_col <- c('#556B2F',"#EE0000","#FF7F24","#8B2500",'#8B1C62',"#436EEE")
  
  df <- df %>% mutate_at(vars(Species_Name),~gsub('_',' ',.))
  
  spec_labs <- c(expression(bold('Acetaldehyde')),expression(bold('Acetic Acid')),expression(bold('Acetic Acid')),
                 expression(bold('Acetonitrile')),expression(bold('Benzene')),expression(bold(C[10]~Aromatics)),
                 expression(bold(C[8]~Aromatics)),expression(bold(C[9]~Aromatics)),expression(bold(CH[4])),
                 expression(bold(Cl^{' -'}~Aerosol)),expression(bold('CO')),expression(bold(CO[2])),
                 expression(bold(D[5]~Siloxane)),expression(bold('Ethanol')),
                 expression(bold('Formaldehyde')),expression(bold('Formic Acid')),
                 expression(bold('Furan')),expression(bold('Furfural')),expression(bold('Hexanone')),
                 expression(bold('Hydroxyacetone')),expression(bold('Isoprene')),
                 expression(bold('Maleic Anhydride')),expression(bold('MEK')),expression(bold('Methanol')),
                 expression(bold('Methylfuran')),expression(bold('Methylfurfural')),expression(bold('Monoterpenes')),
                 expression(bold('MVK/MACR')),expression(bold({NH[4]^{' +'}})),expression(bold(NO)),
                 expression(bold(NO[2])),expression(bold({NO[3]^{'-'}~Aerosol})),
                 expression(bold(NO[x])),expression(bold('OA')),expression(bold('Ozone')),
                 expression(bold('Pentanone')),expression(bold({SO[' 4']^{' 2-'}~Aerosol})),
                 expression(bold('TNRM - Aerosol')),expression(bold('Toluene')))
  
  fin_fig_S6 <- ggplot(df,aes(x=Species_Name,y=perc_cont))+
    geom_col(aes(fill = PMF_factor),color = 'black')+
    scale_y_continuous(breaks = seq(0,100,10),
                       ,minor_breaks = seq(0,100,5),expand = c(0,0))+
    scale_fill_manual(values = factor_col,labels = PMF_factors)+
    scale_x_discrete(labels = spec_labs)+
    xlab('')+
    ylab('% of Species Sum')+
    theme_bw()+
    theme(text = element_text(face = 'bold',family = 'serif',color = 'black'),
          legend.position = 'top',legend.title = element_blank(),
          legend.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 12),
          axis.text = element_text(family = 'serif',face = 'bold',color = 'black',size = 10),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5,size = 12),
          axis.title = element_text(family = 'serif',face = 'bold',color = 'black',size = 15))+
    guides(fill = guide_legend(nrow=1))
  
  return(fin_fig_S6)
}
