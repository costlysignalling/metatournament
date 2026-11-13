# Run from the repository root. Each stage runs in a separate R session.
# Run with no argument to attempt all stages or supply one or more stage names.
stages <- list(
  foundations = c('Figure01.R','Figure02.R','FigureS01.R','FigureS02.R',
                 'Figure_Etude2_three_tournaments_ternary.R',
                 '04_known_talents_strategy_map.R','04_known_talents_utility_contours_equal_sums.R',
                 '04_known_talents_utility_contours_equal_sums_diminishing.R',
                 '06_parameter_concordance_guide.R'),
  participation = c("source('reproduce_participation_figures.R'); reproduce_participation_figures(recompute=TRUE)"),
  asymmetric = c('05_many_tournaments_and_asymmetric_linear.R',
    '05_draw_asymmetric_examples.R','05_asymmetric_tournaments_figures_S14_S15.R',
    "source('05_explore_asymmetric_tournaments.R'); explore_asymmetric(.8,.4,1,1,legend=TRUE,png=TRUE,pdf=TRUE,file='Figure_asymmetric_tournaments_explorer_08_04.png')"),
  coexistence = c(
    "source('06_prepare_coexistence_deterministic.R'); prepare_coexistence('beta'); prepare_coexistence('phi')",
    "source('06_refine_coexistence_curves.R'); refine_coexistence_curves('beta'); refine_coexistence_curves('phi')",
    '06_strategy_coexistence_parameter_space_fixed_and_responsive_examples.R',
    '06_strategy_coexistence_phi_space_by_beta_fixed_and_responsive.R',
    "source('05_correct_energy_phispace_population.R'); correct_energy_phispace_population()",
    '05_energy_allocation_phispace.R'),
  equal_talent = c('05_prepare_equal_talent_geometry_inputs.R','05_refine_equal_talent_nash.R',
    '05_refine_equal_talent_conditional.R','05_equal_talent_5050_parameter_space.R','05_equal_talent_5050_phispace.R'),
  energy = c('05_prepare_energy_allocation_inputs.R','05_refine_energy_allocation_boundaries.R',
    '05_energy_allocation_parameter_space.R','05_prepare_energy_allocation_phispace_inputs.R',
    '05_refine_energy_allocation_phispace_boundaries.R'),
  many = c('06B_verify_many_tournaments.R','06B_refine_many_population.R',
    '06B_many_tournament_figures.R','06B_draw_actor_example.R','06B_prepare_cube_web.R'),
  optout = c('07_optout_figures.R','07_paper_optout_examples.R',
    "source('07_prepare_optout_energy.R'); prepare_optout_energy(eta_values=c(.3,.5,.7))",
    "source('07_refine_optout_energy_boundaries.R'); for(e in c(.3,.5,.7))refine_optout_energy_boundaries(eta=e)",
    vapply(c(.3,.5,.7),function(e)sprintf("options(metatournament.eta=%s); source('07_optout_energy_parameter_space.R')",e),"")),
  boredom = c('08_prepare_boredom_figures.R','08_prepare_paper_boredom.R',
    "source('08_boredom_examples.R'); plot_boredom_examples(variant='more_boredom')",
    '08_boredom_energy_parameter_space.R',"source('08_refine_boredom_level.R'); refine_boredom_level()",'08_boredom_level_parameter_space.R'),
  harvesting = c("source('09_reproduce_harvesting_figures.R'); reproduce_harvesting_figures(recalculate=TRUE)",
    '09_harvesting_explanations.R','09_transformed_talent_density.R',
    "source('09_training_comparison.R'); draw_training_comparison()"),
  convincing_examples = c('10_reproduce_convincing.R',
    "source('10_compare_convincing_models.R'); for(s in comparison_sets)comparison_reproduce_set(s)",
    "source('10_comparison_render.R'); for(s in c('symmetric_beta04','symmetric_beta08','asymmetric'))comparison_render_set(s,include_members=FALSE)",
    '10_refine_convincing_example_maps.R',
    "source('10_local_global_prepare_figures.R'); dir.create('10_local_global_results',showWarnings=FALSE); for(s in c('symmetric_beta04','symmetric_beta08','asymmetric'))saveRDS(lapply(readRDS(paste0('10_convincing_results/',s,'.rds')),campaign_entry_curves),paste0('10_local_global_results/entry_local_',s,'.rds')); source('10_refine_convincing_effort_geometry.R'); prepare_campaign_effort_geometries()",
    '10_local_global_prepare_figures.R','10_paper_convincing_figure.R'),
  convincing_grid = c("source('10_compute_convincing_efficiencies.R'); compute_convincing_efficiencies(export_figures=TRUE,workers=as.integer(Sys.getenv('METATOURNAMENT_CAMPAIGN_WORKERS','1')))",
    '10_verify_campaign_precision_maps.R'),
  manuscripts = 'build_manuscripts.R'
)

args <- commandArgs(trailingOnly=TRUE)
full_run <- !length(args)
if(full_run) args <- c('foundations','participation','asymmetric','energy',
 'coexistence','equal_talent','many','optout','boredom','harvesting',
 'convincing_examples','convincing_grid')
stopifnot(all(args %in% names(stages)))
dir.create('reproduction_logs',showWarnings=FALSE)
rscript <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
failed <- character()
for(stage in args){
  for(i in seq_along(stages[[stage]])){
    command <- stages[[stage]][i]
    log <- file.path('reproduction_logs',sprintf('%s-%02d.log',stage,i))
    message(stage,' ',i,': ',command)
    invocation <- if(file.exists(command)) shQuote(command) else c('-e',shQuote(command))
    status <- system2(rscript,c('--vanilla',invocation),stdout=log,stderr=log)
    if(status!=0L){
      warning('Stage ',stage,' failed. See ',log,immediate.=TRUE)
      failed <- c(failed,stage)
      break
    }
  }
}
if(full_run || 'manuscripts' %in% args){
  missing <- Filter(function(path)!file.exists(path)||file.info(path)$size==0,readLines('figures.txt'))
  if(length(missing)){message('Missing figures: ',paste(missing,collapse=', ')); failed<-c(failed,'figure coverage')}
}
if(length(failed))stop('Failed stages: ',paste(failed,collapse=', '))
