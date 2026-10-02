source("global.R")

ui <- fluidPage(

  titlePanel("MINA-TB Trial Design Explorer"),

  sidebarLayout(

    sidebarPanel(
      width = 3,

      h4("Model parameters"),
      selectInput("country", "Population age distribution",
                  choices=country_choices, selected="South Africa"),
      fluidRow(
        column(6, numericInput("minage", "Recruitment min age", value=15, min=0, max=98, step=1)),
        column(6, numericInput("maxage", "Recruitment max age", value=49, min=1, max=99, step=1))
      ),
      selectInput("foi_structure", "Force of infection",
                  choices=foi_choices, selected=foi_choices[1]),
      selectInput("progression_structure", "Probability of being a fast progressor",
                  choices=pslow_choices, selected="Vynnycky-high"),
      selectInput("slow_structure", "Annual rate of progression for slow progressors",
                  choices=muslow_choices, selected="0.003"),

      hr(),
      h4("Simulation"),
      numericInput("reps", "Simulation reps", value=5, min=1, max=20, step=1),

      hr(),
      actionButton("run", "Run simulation", class="btn-primary", width="100%"),
      br(), br(),
      helpText("Approx run time: ~1 second per rep (about 5 seconds for 5 reps).")
    ),

    mainPanel(
      width = 9,

      plotOutput("main_plot", height="350px"),

      hr(),
      fluidRow(
        column(5,
          h4("Trial size"),
          tableOutput("trial_size")
        ),
        column(7,
          h4("Epidemiological benchmarks"),
          tableOutput("benchmarks")
        )
      )
    )
  )
)


server <- function(input, output, session) {

  sim_results <- reactiveVal(NULL)
  sim_results_pars <- reactiveVal(NULL)

  observeEvent(input$run, {

    pars_base <- list(
      minage       = input$minage,
      maxage       = input$maxage,
      rho          = foi_lookup[[input$foi_structure]],
      p_slow       = extract_probability_of_being_slow(progressor_type_by_age_studies, input$progression_structure),
      mu_slow      = define_mu_slow(input$slow_structure),
      mu_fast      = 1.5,
      p_revert     = default_p_revert,
      mu_revert    = default_mu_revert,
      sigma        = 100,
      trial_length = 3
    )

    agedist  <- extract_age_distribution(unwpp, input$country)
    sigmavec <- c(2, 80)
    reps     <- input$reps
    grid     <- expand.grid(sigma=sigmavec, rep=1:reps)
    n_total  <- nrow(grid)
    casetarget <- 50

    results_list <- vector("list", n_total)

    withProgress(message="Running simulations...", value=0, {
      for (i in 1:n_total) {
        pars_i       <- pars_base
        pars_i$sigma <- grid$sigma[i]
        out          <- sim_stoch(pars_i, casetarget=casetarget, agedist=agedist)
        results_list[[i]] <- tibble(
          sigma           = grid$sigma[i],
          rep             = grid$rep[i],
          n_tested        = out$n_tested,
          n_enrolled      = out$n_enrolled,
          n_fast          = out$n_fast,
          n_slow          = out$n_slow,
          n_overall       = out$n_overall,
          n_trial_cases   = out$n_trial_cases,
          inf_prev_15to24 = out$inf_prev_15to24,
          inf_prev_25to34 = out$inf_prev_25to34
        )
        incProgress(1/n_total,
                    detail=paste0("sigma=", grid$sigma[i], ", rep ", grid$rep[i], "/", reps))
      }
    })

    sim_results(bind_rows(results_list))
    # pars_base is fixed per run, needed again (unchanged) by the benchmarks table below
    sim_results_pars(pars_base)
  })


  # Main plot: box + whiskers, enrolled + tested on same axes
  output$main_plot <- renderPlot({
    req(sim_results())

    plot_df <- sim_results() %>%
      mutate(sigma_label = factor(
        ifelse(sigma==2, "TASA (σ=2)", "IGRA (σ=80)"),
        levels=c("TASA (σ=2)", "IGRA (σ=80)")
      )) %>%
      select(sigma_label, rep, n_enrolled, n_tested) %>%
      pivot_longer(cols=c(n_enrolled, n_tested), names_to="metric", values_to="value") %>%
      mutate(metric = case_when(metric=="n_enrolled"~"Enrolled", metric=="n_tested"~"Tested"))

    ggplot(plot_df, aes(x=sigma_label, y=value, fill=metric)) +
      geom_boxplot(position=position_dodge(width=0.35), width=0.3, alpha=0.8, coef=Inf) +
      scale_fill_brewer(palette="Set1") +
      theme_classic() +
      theme(legend.title=element_blank(), strip.background=element_blank()) +
      labs(x=NULL, y="Number (for 50 disease-endpoint cases)")
  })


  # Trial size table: TASA first, 0 decimal places
  output$trial_size <- renderTable({
    req(sim_results())

    sim_results() %>%
      mutate(sigma_label = factor(
        ifelse(sigma==2, "TASA (σ=2)", "IGRA (σ=80)"),
        levels=c("TASA (σ=2)", "IGRA (σ=80)")
      )) %>%
      group_by(sigma_label) %>%
      summarise(
        Tested    = round(median(n_tested)),
        Enrolled  = round(median(n_enrolled)),
        Fast      = round(median(n_fast)),
        Slow      = round(median(n_slow)),
        Overall   = round(median(n_overall)),
        .groups   = "drop"
      ) %>%
      arrange(sigma_label) %>%
      rename(` ` = sigma_label)
  }, digits=0, striped=TRUE, hover=TRUE, bordered=TRUE,
     caption="All values are medians across simulation reps, for 50 disease-endpoint cases.")


  # Epidemiological benchmarks table - reuses the same sim_results() stochastic output (sigma=80,
  # i.e. lifelong/IGRA-like positivity) rather than re-simulating, via the same
  # summarise_case_incidence_and_ARTI() helper used throughout run_analysis.qmd.
  output$benchmarks <- renderTable({
    req(sim_results())
    req(sim_results_pars())

    pars <- sim_results_pars()
    reps_i <- length(unique(sim_results()$rep))

    metrics <- summarise_case_incidence_and_ARTI(sim_results(), pars=pars, reps=reps_i, sig=80)

    tibble(
      Metric = c(
        "Infection prevalence, age 15-24",
        "Infection prevalence, age 25-34",
        "ARTI (median)",
        "Cases per 100k per year (median)"
      ),
      Value = c(
        sprintf("%.2f", median(metrics$inf_prev_15to24)),
        sprintf("%.2f", median(metrics$inf_prev_25to34)),
        sprintf("%.3f", median(metrics$ARTI)),
        sprintf("%.0f", median(metrics$cases))
      )
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

}


shinyApp(ui, server)
