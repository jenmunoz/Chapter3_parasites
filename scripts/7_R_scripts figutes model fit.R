pdf(file="figures/figures_pdf_manuscript/FigureS3a.Infection_sociality_fit.pdf", width =10, height =10)
pp_check(selected_ecto_infection_brms_bayes_no_int, type = "dens_overlay", ndraws = 200) 
dev.off()

plot <- ggplot2::last_plot() + 
theme_classic(30)
ggsave("infection_1.pdf", plot, width = 10, height = 10)



# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3a2.Infection_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_ecto_p_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 200)+ xlim(0, 2)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("infection_2.pdf", plot, width = 10, height = 10)



# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3a3.Infection_sociality_strength_fit.pdf", width =10, height =10)
pp_check(selected_ecto_p_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 200)+ xlim(0,2)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("infection_3.pdf", plot, width = 10, height = 10)

color_scheme_set("teal")


# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3b1.Lice_abundance_sociality_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_lice_brms_bayes_no_int_priors, type = "dens_overlay", ndraws = 200)+ xlim(0, 20)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("abundance_lice_1.pdf", plot, width = 10, height = 10)


# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3b2.Lice_abundance_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_lice_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 200)+ xlim(0, 50)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("abundance_lice_2.pdf", plot, width = 10, height = 10)

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3b3.Lice_abundance_sociality_strength_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_lice_brms_bayes_no_int_strength_prior, type = "dens_overlay", ndraws = 200)+ xlim(0, 50)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classics(30)
ggsave("abundance_lice_3.pdf", plot, width = 10, height = 10)


color_scheme_set("green")
# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3c1.Mites_abundance_sociality_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_prior, type = "dens_overlay", ndraws = 200)+ xlim(0, 20)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("abundance_mites_1.pdf", plot, width = 10, height = 10)



# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3c2.Mites_abundance_sociality_degree_fit.pdf", width =15, height =10)
pp_check(selected_zinb_a_nf_mites_brms_bayes_no_int_degree_prior, type = "dens_overlay", ndraws = 200)+ xlim(0, 10)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("abundance_mites_2.pdf", plot, width = 10, height = 10)

# model fit
pp_check(ecto_p_brms_bayes_no_int_species_priors_zobi, type = "dens_overlay", ndraws = 200)+ xlim(0, 5)
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("abundance_mites_3.pdf", plot, width = 10, height = 10)

color_scheme_set("orange")
# model fit
pp_check(ecto_p_brms_bayes_no_int_species_priors_degree_zobi, type = "dens_overlay", ndraws = 200)+ xlim(0, 5)
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("prevalence1.pdf", plot, width = 10, height = 10)

# model fit
pp_check(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi, type = "dens_overlay", ndraws = 200)+ xlim(0, 5)
plot <- ggplot2::last_plot() + 
  theme_classic(30)
ggsave("prevalence2.pdf", plot, width = 10, height = 10)

pdf(file="figures/figures_pdf_manuscript/FigureS3d2.Prevalence_socialitystrength_fit.pdf", width =15, height =10)
pp_check(ecto_p_brms_bayes_no_int_species_priors_w_degree_zobi, type = "dens_overlay", ndraws = 200)+ xlim(0, 5)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classics(30)
ggsave("prevalence3.pdf", plot, width = 10, height = 10)

color_scheme_set("red")

# model fit
pdf(file="figures/figures_pdf_manuscript/FigurS3e1.LiceRichness_sociality_fit.pdf", width =15, height =10)
pp_check(selected_poisson_lice_diversity_sociality_no_int_priors, type = "dens_overlay", ndraws = 200)+ xlim(0, 20)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classics(30)
ggsave("rich1.pdf", plot, width = 10, height = 10)

pdf(file="figures/figures_pdf_manuscript/FigureS3e2.LiceRichness_socialitydegree_fit.pdf", width =15, height =10)
pp_check(selected_poisson_lice_diversity_degree_no_int_priors, type = "dens_overlay", ndraws = 200)+ xlim(0, 20)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classics(30)
ggsave("rich2.pdf", plot, width = 10, height = 10)

# model fit
pdf(file="figures/figures_pdf_manuscript/FigureS3e3.LiceRichness_socialityStrength_fit.pdf", width =15, height =10)
pp_check(selected_poisson_lice_diversity_w_degree_no_int_priors, type = "dens_overlay", ndraws = 200)+ xlim(0, 20)
dev.off()
plot <- ggplot2::last_plot() + 
  theme_classics(30)
ggsave("rich3.pdf", plot, width = 10, height = 10)
