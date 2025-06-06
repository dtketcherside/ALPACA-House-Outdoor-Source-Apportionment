# Figure 1: Time Series --------

figure_1 <- gen_figure_1(df=hr_data,CTC_df = CTC_data)

# Figure 2: O3/NO3 Production Windrose --------

# # Only Needs to Be Run if Recreating Figure 2b prior to cropping whitespace
# figure_2b <- gen_figure_2b(df=CTC_data,k_val = 100)
# 
# png(filename = 'output/Figure 2b_alt.png',width = 5,height = 7,res = 1200,units = 'in')
# print(figure_2b)
# dev.off()

figure_2a <- gen_figure_2a(df=CTC_data)

fig2b_gg <- img2gg('output/Figure 2b.png')

figure_2c <- gen_figure_2c(df=NO3_prod_data)

figure_2d <- gen_figure_2d(df=NO3_rad_data)

fig2_final <- wrap_elements(
  ((figure_2a + plot_spacer() + fig2b_gg + plot_layout(widths = c(1,0.15, 1))) / 
     (figure_2c + figure_2d) & 
     theme(plot.tag.position  = c(0.2, .96))) +
    plot_annotation(tag_levels = "a", tag_suffix = ')') +
    plot_layout(heights = c(1.25, 1))
) +
  theme(plot.tag = element_text(size = rel(1), face = 'bold', family = 'serif'))

# Figure 3: Diel Profiles of Key Species --------

figure_3 <- gen_figure_3(df=hr_data)

# Figure 4: Scatterplots --------

figure_4 <- gen_figure_4(df=hr_data)

# Figure 5: Pie Charts of ROC --------

figure_5 <- gen_figure_5(df=hr_data,met_df = met_data)

# Figure 6: PMF Diurnal Profiles and Speciated Contributions --------

figure_6 <- wrap_elements(((gen_figure_6d(df=PMF_TS_data) + theme(plot.tag.position  = c(0.04, .95)) )/
                             ((gen_figure_6b(df=PMF_TS_data) + theme(plot.tag.position  = c(0.08, .95)))+
                              (gen_figure_6a(df=PMF_TS_data) + theme(plot.tag.position  = c(0.1, .95))))/
                             (gen_figure_6c(df=PMF_TS_data) + theme(plot.tag.position  = c(0.04, .95))))+
                            plot_annotation(tag_levels = "a",tag_suffix = ')')+
                            theme(plot.tag = element_text(size = rel(1),face = 'bold',family = 'serif')) +
                            plot_layout(heights = c(5, 10, 5)) & ylab(NULL)) +
  labs(tag = expression(bold("Normalized g-values ("*bar(x)*' = 1)'))) +
  theme(
    plot.tag = element_text(face = 'bold',family = 'serif',size = rel(1.1), angle = 90),
    plot.tag.position = "left"
  )

# Figure 7. PMF Factor Contributions to Species Concentrations ----------

figure_7 <- gen_figure_7(tabs4)

################### Save Main Text Figures #############################

ggsave(filename = 'output/Figure 1.tiff',
       plot = figure_1,width = 1000,height = 900,scale = 3,units = 'px')

ggsave(filename = 'output/Figure 2.tiff',
       plot = fig2_final,width = 1000,height = 750,scale = 3,units = 'px')

ggsave(filename = 'output/Figure 3.tiff',
       plot = figure_3,width = 1000,height = 900,scale = 3,units = 'px')

ggsave(filename = 'output/Figure 4.tiff',
       plot = figure_4,width = 950,height = 700,scale = 3,units = 'px')

ggsave(filename = 'output/Figure 5.tiff',
       plot = figure_5,width = 550,height = 500,scale = 3,units = 'px')

ggsave(filename = 'output/Figure 6.tiff',
       plot = figure_6,width = 800,height = 800,scale = 3,units = 'px')

ggsave(filename = 'output/Figure 7.tiff',
       plot = figure_7,width = 750,height = 500,scale = 3,units = 'px')

################### Supplemental Figures ###################

# Anticorrelations ---------

figure_S1 <- gen_figure_1_SI(df = CTC_data)

ggsave(filename = 'output/Supplemental/Figure S1.tiff',
       plot = figure_S1,width = 1000,height = 500,scale = 3,units = 'px')

# Toluene/Benzene Ratio -----------

figure_S2 <- gen_figure_2_SI(hr_data)

ggsave(filename = 'output/Supplemental/Figure S2.tiff',
       plot = figure_S2,width = 1000,height = 900,scale = 1,units = 'px')


# ADEC Local Source Profile Comparison with PMF Factors -----------

figure_S3 <- gen_figure_3_SI(SI_plot_data)

ggsave(filename = 'output/Supplemental/Figure S3.tiff',
       plot = figure_S3,width = 1200,height = 1600,scale = 2,units = 'px')

# Sulfate Diel Profile by Month --------------

figure_S4 <- gen_figure_4_SI(hr_data)

ggsave(filename = 'output/Supplemental/Figure S4.tiff',
       plot = figure_S4,width = 1500,height = 900,scale = 1,units = 'px')

# PMF Diel Profiles and Time Series ---------

figure_S5 <- gen_figure_5_SI(PMF_TS_data)

ggsave(filename = 'output/Supplemental/Figure S5.tiff',
       plot = figure_S5,width = 1255,height = 783,scale = 3,units = 'px')

# PMF Factor Contributions to Species Sum ---------

figure_S6 <- gen_figure_6_SI(tabs4)

ggsave(filename = 'output/Supplemental/Figure S6.tiff',
       plot = figure_S6,width = 1500,height = 900,scale = 3,units = 'px')
