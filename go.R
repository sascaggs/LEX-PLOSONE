#####################################################
#
# STRAND modeling of reported labor ties 
# Shane A. Scaggs
# PLoS One 
# 2025
#
#####################################################

# Setup ----
# install and load R packages 
#install.packages(c('tidyverse','tidygraph','ggraph','igraph','latex2exp','patchwork','posterior'))
#devtools::install_github('ctross/STRAND')
library(tidyverse) # data wrangling, plotting in ggplot2
library(tidygraph) # network graph wrangling
library(ggraph)    # network graph plotting 
library(igraph)    # create graph object, community mutuality (reciprocity)
library(latex2exp) # latex in ggplot2 titles
library(patchwork) # plot layouts and labels 
library(STRAND)    # latent network modeling 


graph_theme = theme(
     # panels
     panel.background = element_rect(
          color  = 'black', 
          fill   = '#ffffffff', 
          size   = 1 ), 
     panel.grid = element_blank( ), 
     panel.spacing = unit(15, 'pt'), 
     # axes
     axis.ticks  = element_line(
          color = 'black', 
          size  = 0.5 ), 
     axis.ticks.length = unit(2, 'mm'), 
     # strips
     strip.background = element_rect(
          color = '#ffffffff', 
          fill  = '#ffffffff',), 
     strip.text = element_text(
          color  = 'black', 
          vjust  = 1.2, 
          size   = 10, 
          margin = unit( c(2,0,2,0), 'mm') ), 
     # title
     plot.title = element_text(
          size  = 12, 
          hjust = 0.5 )
)

theme_set(graph_theme)

# outcome network: layer 1 = who would you ask for labor, layer 2 = who would you give labor to, if asked
Ask4Labor = as.matrix( read.csv('Data/Ask4Labor.csv', header = T) )
rownames(Ask4Labor) = colnames(Ask4Labor)

GiveLabor = as.matrix( read.csv('Data/GiveLabor.csv', header = T) )
rownames(GiveLabor) = colnames(GiveLabor)

outcome = list(
     Ask4Labor = Ask4Labor, 
     GiveLabor = GiveLabor
)


# dyad covariates
# Relatedness: from kinship2 pedigrees based on validated kin relations
# Church: binary indicator of whether two farmers attend the same church 
Relatedness = as.matrix( read.csv('Data/Relatedness.csv', header = T) )
rownames(Relatedness) = colnames(Relatedness)
SameChurch = as.matrix( read.csv('Data/SameChurch.csv', header = T) )
rownames(SameChurch) = colnames(SameChurch)

dyad = list(
     Relatedness = Relatedness, 
     SameChurch  = SameChurch
)

# individual covariates
# Material wealth and cattle pasture land use 
indiv = read.csv('Data/indiv.csv', header = T)
rownames(indiv) = rownames(Ask4Labor)

# group covariates
# Village of residence 
groups = read.csv('Data/groups.csv', header = T)
groups = data.frame(
     Church  = as.factor(groups$Church), 
     Village = as.factor(groups$Village)
)
rownames(groups) = rownames(Ask4Labor)

# make clean strand data 
dat = make_strand_data(
     outcome = outcome,
     individual_covariates = indiv, 
     dyadic_covariates = dyad,
     block_covariates = groups,
     outcome_mode = "bernoulli", 
     link_mode = "logit", 
     imputation = TRUE
)


# Plot networks ----

## Step 1 ---------------
# layer 1
ask_g = graph_from_adjacency_matrix(outcome$Ask4Labor)

# vertex and edge attributes
V(ask_g)$Village = groups$Village
E(ask_g)$Reciprocal = which_mutual(ask_g)

# layer 2
give_g = graph_from_adjacency_matrix(outcome$GiveLabor)

# vertex and edge attributes
V(give_g)$Village = groups$Village
E(give_g)$Reciprocal = which_mutual(give_g)

## Step 2 ---------------
# layout 
set.seed(666)
nice = layout_nicely(ask_g)

# network graphs 
nets = ask_g |> 
     as_tbl_graph() |> 
     ggraph(layout = nice) + 
     geom_edge_link( aes(color = Reciprocal) ) + 
     geom_node_point( aes(fill = Village, color = Village), 
                      shape  = 21, 
                      stroke = 0.681 ) + 
     theme_graph() + 
     theme( legend.position = 'none', 
            plot.margin = unit(c(0,0,0,0),'cm'), 
            plot.title = element_text(
                 size  = 14, 
                 hjust = 0.5, 
                 face  = 'plain'), 
            plot.subtitle = element_text(
                 size  = 12, 
                 hjust = 0.5, 
                 face  = 'italic') ) + 
     scale_edge_color_manual(
          values = c('#00000035','red') ) + 
     scale_fill_manual(
          values = c('#8999ff',
                     '#44AA99',
                     '#CC6677',
                     '#ffddaa',
                     '#fff') ) +
     scale_color_manual(
          values = c('#332288',
                     '#117733',
                     '#882255',
                     '#DDCC77',
                     '#666') ) +
     ggtitle('Layer 1', 
             subtitle = "Who would you ask for labor?" ) + 
     
     give_g |> 
     as_tbl_graph() |> 
     ggraph(layout = nice) + 
     geom_edge_link( aes(color = Reciprocal) ) + 
     geom_node_point( aes(fill = Village, color = Village), 
                      shape  = 21, 
                      stroke = 0.681 ) + 
     theme_graph() + 
     theme( plot.margin = unit(c(0,0,0,0),'cm'), 
            plot.title = element_text(
                 size  = 14, 
                 hjust = 0.5, 
                 face  = 'plain'), 
            plot.subtitle = element_text(
                 size  = 12, 
                 hjust = 0.5, 
                 face  = 'italic') ) + 
     scale_edge_color_manual(
          values = c('#00000035','red') ) + 
     scale_fill_manual(
          values = c('#8999ff',
                     '#44AA99',
                     '#CC6677',
                     '#ffddaa',
                     '#fff') ) +
     scale_color_manual(
          values = c('#332288',
                     '#117733',
                     '#882255',
                     '#DDCC77',
                     '#666') ) +
     ggtitle('Layer 2', 
             subtitle = "Who would you give labor to, if asked?")  

## Step 3 -----------------

deg = rbind(
     data.frame(
          layer = 'ask', 
          out_degree = strength(ask_g, mode = 'out'), 
          in_degree  = strength(ask_g, mode = 'in')
     ),
     data.frame(
          layer = 'give', 
          out_degree = strength(give_g, mode = 'out'), 
          in_degree  = strength(give_g, mode = 'in')
     )
) |> 
     gather(key = metric, value = value, -layer) 


## Step 4 -----------------
ask_strip_labs = as_labeller(c(
     `out_degree` = "Out-degree\n(ask for labor)", 
     `in_degree`  = "In-degree\n(ask for labor)")
)

give_strip_labs = as_labeller( c(
     `out_degree` = "Out-degree\n(give labor)", 
     `in_degree`  = "In-degree\n(give labor)")
)

deg_plot = deg |> 
     filter( layer == 'ask' ) |> 
     ggplot( aes(x=value)) + 
     geom_histogram( binwidth = 1, 
                     color = 'black', 
                     fill  = 'white' ) + 
     theme( legend.position = 'none' ) + 
     scale_x_continuous( limits = c(-1,17) ) + 
     facet_wrap(~metric, 
                nrow = 1, 
                labeller = ask_strip_labs) + 
     labs(x=NULL, y='Count') +
     
     deg |> 
     filter( layer == 'give' ) |> 
     ggplot( aes(x=value)) + 
     geom_histogram( binwidth = 1, 
                     color = 'black', 
                     fill  = 'white' ) + 
     theme( legend.position = 'none' ) + 
     scale_x_continuous( limits = c(-1,17) ) + 
     facet_wrap( ~metric, 
                 nrow = 1, 
                 labeller = give_strip_labs ) + 
     labs(x=NULL, y='Count')

## Step 5 -----------------
nets / deg_plot + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') + plot_layout(heights = c(3,1.854))

#ggsave(filename = "Figs/Raw/ObservedNetworkLayers.png", height = 4*1.5, width = 6*1.5, units = 'in', dpi = 1200)


# Fit models ----

## Pasture only ----
run_models = FALSE # set to TRUE to run
if(run_models) {
     fitP = fit_latent_network_model_missings(
          data = dat, 
          block_regression  = ~ Village, 
          focal_regression  = ~ Pasture, 
          target_regression = ~ Pasture, 
          dyad_regression   = ~ Relatedness + SameChurch , 
          rtt_regression    = ~ 1, 
          fpr_regression    = ~ 1, 
          theta_regression  = ~ 1, 
          mode = 'mcmc', 
          return_predicted_network = FALSE, 
          stan_mcmc_parameters = list(
               seed = 666, 
               chains = 4, 
               parallel_chains = 4, 
               refresh = 100, 
               iter_warmup = 1000, 
               adapt_delta = 0.9, 
               iter_sampling = 2000, 
               max_treedepth = NULL )
     )
}

# Save and load output files and STRAND model object

#fitP$fit$save_output_files(dir = 'Fit models')
#save(fitP, file = 'Fit models/fitP.Rdata')
#load("Fit models/fitP.Rdata")

## Pasture + Wealth ----
if(run_models){
     fitPW = fit_latent_network_model_missings(
          data = dat, 
          block_regression  = ~ Village, 
          focal_regression  = ~ Pasture + Wealth, 
          target_regression = ~ Pasture + Wealth, 
          dyad_regression   = ~ Relatedness + SameChurch , 
          rtt_regression    = ~ 1, 
          fpr_regression    = ~ 1, 
          theta_regression  = ~ 1, 
          mode = 'mcmc', 
          return_predicted_network = FALSE, 
          stan_mcmc_parameters = list(
               seed = 666, 
               chains = 4, 
               parallel_chains = 4, 
               refresh = 100, 
               iter_warmup = 1000, 
               adapt_delta = 0.9, 
               iter_sampling = 2000, 
               max_treedepth = NULL)
     )
}

# Save and load output files and STRAND model object

#fitPW$fit$save_output_files(dir = 'Fit models')
#save(fitPW, file = 'Fit models/fitPW.Rdata')
load('Fit models/fitPW.Rdata')

## Diagnostics ---- 

# check divergence, e-bfmi, exceeding max treedepth 
# draws and rhat  
see_diagnostics = FALSE # requires a model fit 
if(see_rhat){
     
     fitP$fit$diagnostic_summary()
     fitPW$fit$diagnostic_summary()
     
     drawsP = fitP$fit$draws(format = 'df')
     summaryP = posterior::summarise_draws(drawsP)
     max(summaryP$rhat, na.rm = T)
     summaryP |> ggplot(aes(x=rhat)) + geom_histogram(bins = 100)
     
     drawsPW = fitPW$fit$draws(format = 'df')
     summaryPW = posterior::summarise_draws(drawsPW)
     max(summaryPW$rhat, na.rm = T)
     summaryPW |> ggplot(aes(x=rhat)) + geom_histogram(bins = 100)
}

# Table 2 Results ----

# These are the STRAND plots shown in ctross/STRAND, useful for quick inspection and looking at measurement model

## P Results ----  
 
resP = summarize_strand_results(fitP) # requires a model fit
resP$summary
visP1 = strand_caterpillar_plot(resP, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE) 
visP1
visP2 = strand_caterpillar_plot(resP, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), normalized=TRUE)
visP2
visP3 = strand_caterpillar_plot(resP, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), only_slopes=FALSE, normalized=FALSE, only_technicals=TRUE)
visP3


## PW Results ----  
resPW = summarize_strand_results(fitPW) # requires a model fit
resPW$summary
visPW1 = strand_caterpillar_plot(resPW, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE) 
visPW1
visPW2 = strand_caterpillar_plot(resPW, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), normalized=TRUE)
visPW2
visPW3 = strand_caterpillar_plot(resPW, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), only_slopes=FALSE, normalized=FALSE, only_technicals=TRUE)
visPW3


## Reciprocity Correlations ----
cor1 = strand_VPCs(fitP, n_partitions = 4, include_reciprocity = T, mode = 'cor') 
cor1$Reciprocity

cor2 = strand_VPCs(fitPW, n_partitions = 4, include_reciprocity = T, mode = 'cor')
cor2$Reciprocity


# Generate Figures using the Pasture+Wealth model ---- 
fit2 = fitPW$fit

## Block effects ----

# get scalar block for baseline 
block_mat = resPW$samples$srm_model_samples$block_parameters[[1]]
block_mean = mean(block_mat)

# get blocks by village
block_matrix = apply(resPW$samples$srm_model_samples$block_parameters[[2]], c(2,3), mean )

# posterior probability and sd
block_matrix_p = apply(resPW$samples$srm_model_samples$block_parameters[[2]], c(2,3), function(x) mean( inv_logit(x)) )
block_matrix_sd = apply(resPW$samples$srm_model_samples$block_parameters[[2]], c(2,3), function(x) sd( inv_logit(x)))

# block effects labels
colnames(block_matrix) = c('CS','GC','LS','MA','U')
rownames(block_matrix) = c('CS','GC','LS','MA','U')
colnames(block_matrix_p) = c('CS','GC','LS','MA','U')
rownames(block_matrix_p) = c('CS','GC','LS','MA','U')
colnames(block_matrix_sd) = c('CS','GC','LS','MA','U')
rownames(block_matrix_sd) = c('CS','GC','LS','MA','U')

# marginalize over directional effects of blocks (e.g., V_i -> V_j vs V_- <- V_j)
sym_mean = (block_matrix_p + t(block_matrix_p)) / 2
sym_sd = (block_matrix_sd + t(block_matrix_sd)) / 2

get_lower_tri <- function(mat) {
     mat[upper.tri(mat)] <- NA
     mat
}

# pairwise village probabilities 
mean_df <- reshape2::melt(get_lower_tri(sym_mean), na.rm = TRUE)
sd_df   <- reshape2::melt(get_lower_tri(sym_sd), na.rm = TRUE)

df <- cbind(mean_df, sd = sd_df$value)
colnames(df) <- c("Vi", "Vj", "mean", "sd")

# plot mean and sd 
df |> 
     ggplot( aes(y = paste(Vi, Vj, sep='-'), x = mean) ) + 
     geom_vline( aes(xintercept = inv_logit(block_mean)), lty=2 )  + 
     geom_errorbar( aes(xmin = mean-sd, xmax = mean+sd), width=0 ) + 
     geom_point( fill='white', pch=21 ) + 
     labs(x = TeX("$\\textit{\\bar{p}_{V_iV_j}} \\ \\pm \\ \\textit{\\sigma_{p}}$"), 
          y = TeX("$\\textit{V_iV_j}"))

#ggsave(filename = "Figs/Raw/BlockEffects.png", height = 3, width = 4, units = 'in', dpi = 1200)

# Figure 4 - Dyadic effects ---- 

# pull out data and imputed 
dimnames(fitPW$data$dyad_set)[[3]]
fitPW$data$locations_missing_dyad_set

# imputed values 
imputed_vals = fit2$draws(variables = 'imp_dyad_set', format = 'draws_matrix')
imp_means = colMeans(imputed_vals)

# original dyadic array
dyad_array_imputed = fitPW$data$dyad_set

# fill in missing values 
locations = fitPW$data$locations_missing_dyad_set
for(k in seq_along(imp_means)) {
     i = locations[k,1]
     j = locations[k,2]
     cov_idx = locations[k,3]
     dyad_array_imputed[i,j, cov_idx] = imp_means[k]
}

# covariate matrices
relatedness_mat = dyad_array_imputed[,,2]
colnames(relatedness_mat) = 1:ncol(relatedness_mat)
church_mat = dyad_array_imputed[,,3]
colnames(church_mat) = 1:ncol(church_mat)

# converting to dyadic format with imputation indicators  
rel_df = as_tibble(relatedness_mat) |> 
     mutate(i = row_number()) |> 
     pivot_longer(-i, names_to = 'j', values_to = 'Relatedness') |> 
     mutate(j = as.integer(j)) |> 
     filter(i != j)

church_df = as_tibble(church_mat) |> 
     mutate(i = row_number()) |> 
     pivot_longer(-i, names_to = 'j', values_to = 'SameChurch') |> 
     mutate(j = as.integer(j)) |> 
     filter(i != j)

imputed_locs_rel = as_tibble(locations) |> 
     filter(dim3 == 2) |> 
     transmute(i = dim1, j = dim2, rel_was_imputed = T) 

imputed_locs_ch = as_tibble(locations) |> 
     filter(dim3 == 3) |> 
     transmute(i = dim1, j = dim2, church_was_imputed = T) 

# final 
dyads_df = merge(rel_df, church_df, by = c('i','j'))

dyads_df = dyads_df |> 
     left_join(imputed_locs_rel, by = c('i','j')) |> 
     mutate(rel_was_imputed = replace_na(rel_was_imputed, FALSE)) |> 
     left_join(imputed_locs_ch, by = c('i','j')) |> 
     mutate(church_was_imputed = replace_na(church_was_imputed, FALSE)) |> 
     arrange(i,j)


# get dyadic random effects and store in tidy format
# dyadic effects 
dyad_re = resPW$samples$srm_model_samples$dyadic_random_effects

# get point estimate
dyad_re_mean = apply(dyad_re, c(2,3), mean)
colnames(dyad_re_mean) = 1:ncol(dyad_re_mean)

# create the dyadic format and join for each direction 
dyad_re_df = as_tibble(dyad_re_mean) |> 
     mutate(i = row_number()) |> 
     pivot_longer(
          -i, names_to = 'j', values_to = 'effect_ij'
     ) |> 
     mutate(j = as.integer(j)) |> 
     filter(i != j) 


dyad_re_df = dyad_re_df |>   
     inner_join(
          dyad_re_df %>% rename(i_rev = j, j_rev = i, effect_ji = effect_ij), by = c("i" = 'i_rev', "j" = 'j_rev')
     ) |> 
     filter( i < j )

# final merge for dyad effects
dyads_df = merge(dyads_df, dyad_re_df, by = c('i','j'))


# add blocks as columns and merge 
# to do this we need the generalized indices since these hold the village membership information 
gen_array = resPW$samples$srm_model_samples$focal_target_random_effects

gen_mean = apply(gen_array, c(2,3), mean)

gen_df = tibble(
     id = 1:155,
     give = gen_mean[,1], 
     receive = gen_mean[,2]
)

# get the block intercepts and attach to each individual i 
block_ids = fitPW$data$block_set[,2]
block_mean_by_i = sapply(1:155, function(i) { mean(block_matrix[block_ids[i], ]) })
block_df = data.frame(i = 1:155, block_ids)
gen_df$block_id = block_df$block_ids

# merge to get all giving and receiving effects 
dy_tmp = merge(dyads_df, gen_df, by.x = 'i', by.y = 'id')
colnames(dy_tmp)[9:11] = c('give_eff_i','receive_eff_i','block_id_i')
dyads_df_block = merge(dy_tmp, gen_df, by.x = 'j', by.y = 'id')
colnames(dyads_df_block)[12:14] = c('give_eff_j','receive_eff_j','block_id_j')

# get the block intercept estimates
dyads_df_block$block_effect = mapply(function(bi,bj) block_matrix[bi,bj], 
                                     dyads_df_block$block_id_i, 
                                     dyads_df_block$block_id_j)

# save coef estimates for relatedness and church
dyad_coefs = resPW$samples$srm_model_samples$dyadic_coeffs


# Show marginal effect of relatedness with imputed values  
dyads_df_block |> 
     mutate(
          village_block   = ifelse(block_id_i == block_id_j, 1, 0), 
          SameChurch_draw = rbinom(nrow(dyads_df_block), 1, prob = SameChurch), 
          Kin             = ifelse(Relatedness < 0, 'Non-kin', 'Kin'), 
          # constructed dyadic effects, no church effect  
          prob_ij = inv_logit(
               block_effect + effect_ij + Relatedness * dyad_coefs[1] + SameChurch * 0 ), 
          prob_ji = inv_logit(
               block_effect + effect_ji + Relatedness * dyad_coefs[1] + SameChurch * 0 )) |> 
     ggplot( aes(x = prob_ij, y = prob_ji) ) + 
     geom_point( aes(color = Kin, shape = rel_was_imputed) ) + 
     geom_abline(linetype = 2) + 
     scale_shape_manual(values = c(1,4)) + 
     scale_color_manual(values = c('magenta3','black'), 
                        labels = c(TeX('$r > 0$'), TeX('$r = 0$'))) +  
     theme(legend.position = c(0.2,0.8), 
           legend.key      = element_blank()) + 
     labs(x     = TeX("Dyadic effect ($D_{\\[i,j\\]}$)"), 
          y     = TeX("Dyadic effect ($D_{\\[j,i\\]}$)"), 
          color = 'Dyadic\nrelatedness', 
          shape = 'Relatedness\nimputed') 


# Show marginal effect of church with imputed values  
dyads_df_block |> 
     mutate(
          village_block   = ifelse(block_id_i == block_id_j, 1, 0), 
          SameChurch_draw = rbinom(nrow(dyads_df_block), 1, prob = SameChurch), 
          Kin             = ifelse(Relatedness < 0, 'Non-kin', 'Kin'), 
          # constructed dyadic effects, no church effect  
          prob_ij = inv_logit(
               block_effect + effect_ij + Relatedness * 0 + SameChurch * dyad_coefs[2] ), 
          prob_ji = inv_logit(
               block_effect + effect_ji + Relatedness * 0 + SameChurch * dyad_coefs[2] )) |> 
     ggplot( aes(x = prob_ij, y = prob_ji) ) + 
     geom_point( aes(color = factor(SameChurch_draw), shape = church_was_imputed) ) + 
     geom_abline(linetype = 2) + 
     scale_shape_manual(values = c(1,4)) + 
     scale_color_manual(values = c('black','magenta3'), 
                        labels = c('No','Yes')) +  
     theme(legend.position = c(0.2,0.8), 
           legend.key      = element_blank()) + 
     labs(x     = TeX("Dyadic effect ($D_{\\[i,j\\]}$)"), 
          y     = TeX("Dyadic effect ($D_{\\[j,i\\]}$)"), 
          color = 'Same church', 
          shape = 'Church co-membership\nimputed') 

# plot both marginal effects  
two_panel = dyads_df_block |> 
     mutate(
          village_block   = ifelse(block_id_i == block_id_j, 1, 0), 
          SameChurch_draw = rbinom(nrow(dyads_df_block), 1, prob = SameChurch), 
          Kin             = ifelse(Relatedness < 0, 'Non-kin', 'Kin'), 
          # no church effect
          prob_ij = inv_logit(block_effect + effect_ij + Relatedness * dyad_coefs[1] + 0 * dyad_coefs[2]), 
          prob_ji = inv_logit(block_effect + effect_ji + Relatedness * dyad_coefs[1]+ 0 * dyad_coefs[2])) |> 
     
     ggplot( aes(x = prob_ij, y = prob_ji) ) + 
     geom_point( aes(color = Kin), pch=21 ) + 
     geom_abline( linetype = 2 ) + 
     scale_color_manual( values = c('magenta2','black') ) +  
     xlim( c(0,1) ) + 
     ylim( c(0,1) ) + 
     theme( legend.position = c(0.25,0.8), 
            legend.key      = element_blank()) + 
     # axis titles edited in photoshop 
     labs( x     = TeX("Dyadic effect ($D_{\\[i,j\\]}$)"), 
           y     = TeX("Dyadic effect ($D_{\\[j,i\\]}$)"), 
           fill = 'Kinship' ) +
     
     dyads_df_block |> 
     mutate(
          village_block   = ifelse(block_id_i == block_id_j, 1, 0), 
          SameChurch_draw = rbinom(nrow(dyads_df_block), 1, prob = SameChurch), 
          Church          = ifelse(SameChurch_draw == 0, 'Different church', 'Same church'), 
          Kin             = ifelse(Relatedness < 0, 'Non-kin', 'Kin'), 
          
          # no kin effect
          prob_ij = inv_logit(block_effect + effect_ij + 0 * dyad_coefs[1] + SameChurch * dyad_coefs[2]), 
          prob_ji = inv_logit(block_effect + effect_ji + 0 * dyad_coefs[1]+ SameChurch * dyad_coefs[2])) |> 
     
     ggplot( aes(x = prob_ij, y = prob_ji) ) + 
     geom_point( aes(color = Church), pch=21 ) + 
     geom_abline(linetype = 2) + 
     scale_color_manual(values = c('black','magenta3')) +  
     theme( legend.position = c(0.325,0.8), 
            legend.key      = element_blank() ) + 
     xlim( c(0,1) ) + 
     ylim( c(0,1) ) + 
     # axis titles edited in photoshop 
     labs( x    = TeX("Dyadic effect $\\ \\textit{i \\rightarrow j}$"), 
           y    = TeX("Dyadic effect $\\ \\textit{j \\rightarrow i}$"), 
           fill = 'Church membership' )  

# arrange panels, axis titles are added in photoshop due to LaTeX mathtype limitations in ggplot2
two_panel + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

#ggsave(filename = "Figs/Raw/DyadicEffects.png", height = 3*1.3, width = 6*1.3, units = 'in', dpi = 1200)

## Generalized effects ---- 

# add any-to-any block as baseline 
# add village specific block intercepts 
gen_df$baseline = block_mean
gen_df$village_block = block_mean_by_i

# get imputed pasture values 
focal_imputed_vals = fit2$draws(variables = 'imp_focal_set', format = 'draws_matrix')
focal_imputed_mean = colMeans(focal_imputed_vals)

gen_df$Pasture = fitPW$data$focal_set[,2]
gen_df[ which(gen_df$Pasture == -9999999), 'Pasture'] = focal_imputed_mean
gen_df$Wealth = fitPW$data$focal_set[,3]

# pasture effects 
focal_coeffs  = apply(resPW$samples$srm_model_samples$focal_coeffs, 2, mean)
target_coeffs = apply(resPW$samples$srm_model_samples$target_coeffs, 2, mean)

# show imputed values 
gen_df |> 
     mutate(prob_give    = inv_logit(baseline + give + Pasture * focal_coeffs[1]),  
            prob_receive = inv_logit(baseline + receive + Pasture * target_coeffs[1] ), 
            Pasture_imp  = ifelse(Pasture == 0, 'Does not raise cattle', 
                                  ifelse(Pasture == 1, 'Does raise cattle','Imputed'))) |> 
     ggplot( aes( x = prob_give, y = prob_receive) ) + 
     geom_abline( intercept = 0, slope = 1, linetype = 2 ) + 
     geom_point( aes(shape = Pasture_imp) )  + 
     facet_wrap( ~Pasture_imp ) + 
     xlim( c(0,0.006) ) + 
     ylim( c(0,0.006) ) + 
     scale_shape_manual(values = c(21,21,4)) + 
     theme( legend.position = 'none', 
            legend.key = element_blank() ) + 
     labs( x     = TeX('Generalized giving ($\\mathcal{G}_{\\[i\\]}$)'), 
           y     = TeX('Generalized receiving ($\\mathcal{R}_{\\[j\\]}$)'), 
           color = '', 
           fill  = '', 
           shape = '' ) 

# associated between giving and receiving hold wealth constant at the average (0)
# axis titles are added in photoshop due to limitations of LaTeX mathtype in ggplot2
gen_df |> 
     mutate(prob_give    = inv_logit(baseline + give + Pasture * focal_coeffs[1]),  
            prob_receive = inv_logit(baseline + receive + Pasture * target_coeffs[1] ), 
            Pasture_imp  = ifelse(Pasture == 0, 'Does not raise cattle',
                                  ifelse(Pasture == 1, 'Does raise cattle','Imputed'))) |> 
     ggplot( aes( x = prob_give, 
                  y = prob_receive) ) + 
     geom_abline( intercept = 0, slope = 1, linetype = 2 ) + 
     geom_point( aes(color = Pasture_imp, fill = Pasture_imp), shape = 21 )  + 
     xlim( c(0,0.006) ) + 
     ylim( c(0,0.006) ) + 
     scale_color_manual( values = c('black','red','black') ) + 
     scale_fill_manual(  values = c('black','red','#ffffff') ) + 
     theme( legend.position   = c(0.25,0.85), 
            legend.key        = element_blank(), 
            legend.background = element_blank()) +
     labs( x     = TeX('Generalized giving ($\\mathcal{G}_{\\[i\\]}$)'), 
           y     = TeX('Generalized receiving ($\\mathcal{R}_{\\[i\\]}$)'), 
           color = '', 
           fill  = '' ) 

#ggsave(filename = "Figs/Raw/GeneralizedEffects.png", height = 4, width = 4.5, units = 'in', dpi = 1200)

  