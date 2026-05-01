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
        column(6, numericInput("minage", "Recruitment min age", value=18, min=0, max=98, step=1)),
        column(6, numericInput("maxage", "Recruitment max age", value=49, min=1, max=99, step=1))
      ),
      selectInput("foi_structure", "Force of infection",
                  choices=foi_choices, selected="base"),
      selectInput("progression_structure", "Probability of being a fast progressor",
                  choices=pslow_choices, selected="Vynnycky-1997"),
      selectInput("slow_structure", "Annual rate of progression for slow progressors",
                  choices=muslow_choices, selected="base"),

      hr(),
      h4("Simulation"),
      numericInput("reps", "Simulation reps", value=5, min=1, max=20, step=1),

      hr(),
      actionButton("run", "Run simulation", class="btn-primary", width="100%"),
      br(), br(),
      helpText("Approx run time: 20–30 seconds for 5 reps.")
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

  observeEvent(input$run, {

    pars_base <- list(
      minage       = input$minage,
      maxage       = input$maxage,
      rho          = define_foi_by_age(input$foi_structure),
      p_slow       = extract_probability_of_being_slow(progressor_type_by_age_studies, input$progression_structure),
      mu_slow      = define_mu_slow(input$slow_structure),
      mu_fast      = 1.5,
      sigma        = 100,
      trial_length = 3,
      ve           = 0.55
    )

    agedist  <- extract_age_distribution(unwpp, input$country)
    sigmavec <- c(2, 80)
    reps     <- input$reps
    grid     <- expand.grid(sigma=sigmavec, rep=1:reps)
    n_total  <- nrow(grid)

    results_list <- vector("list", n_total)

    withProgress(message="Running simulations...", value=0, {
      for (i in 1:n_total) {
        pars_i       <- pars_base
        pars_i$sigma <- grid$sigma[i]
        out          <- sim_stoch(pars_i, fasttarget=50, agedist=agedist)
        results_list[[i]] <- tibble(
          sigma       = grid$sigma[i],
          rep         = grid$rep[i],
          n_tested    = out$n_tested,
          n_recruited = out$n_recruited,
          n_fast      = out$n_fast,
          n_slow      = out$n_slow,
          n_overall   = out$n_overall
        )
        incProgress(1/n_total,
                    detail=paste0("sigma=", grid$sigma[i], ", rep ", grid$rep[i], "/", reps))
      }
    })

    sim_results(bind_rows(results_list))
  })


  # Main plot: box + whiskers, recruited + screened on same axes
  output$main_plot <- renderPlot({
    req(sim_results())

    plot_df <- sim_results() %>%
      mutate(sigma_label = factor(
        ifelse(sigma==2, "TASA (σ=2)", "IGRA (σ=80)"),
        levels=c("TASA (σ=2)", "IGRA (σ=80)")
      )) %>%
      select(sigma_label, rep, n_recruited, n_tested) %>%
      pivot_longer(cols=c(n_recruited, n_tested), names_to="metric", values_to="value") %>%
      mutate(metric = case_when(metric=="n_recruited"~"Recruited", metric=="n_tested"~"Screened"))

    ggplot(plot_df, aes(x=sigma_label, y=value, fill=metric)) +
      geom_boxplot(position=position_dodge(width=0.35), width=0.3, alpha=0.8, coef=Inf) +
      scale_fill_brewer(palette="Set1") +
      theme_classic() +
      theme(legend.title=element_blank(), strip.background=element_blank()) +
      labs(x=NULL, y="Number (per 50 fast progressors)")
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
        Screened  = round(median(n_tested)),
        Recruited = round(median(n_recruited)),
        Fast      = round(median(n_fast)),
        Slow      = round(median(n_slow)),
        Overall   = round(median(n_overall)),
        .groups   = "drop"
      ) %>%
      arrange(sigma_label) %>%
      rename(` ` = sigma_label)
  }, digits=0, striped=TRUE, hover=TRUE, bordered=TRUE,
     caption="All values are medians across simulation reps.")


  # Epidemiological benchmarks table
  output$benchmarks <- renderTable({
    req(sim_results())

    pars <- list(
      minage       = input$minage,
      maxage       = input$maxage,
      rho          = define_foi_by_age(input$foi_structure),
      p_slow       = extract_probability_of_being_slow(progressor_type_by_age_studies, input$progression_structure),
      mu_slow      = define_mu_slow(input$slow_structure),
      mu_fast      = 1.5,
      sigma        = 80,
      trial_length = 3,
      ve           = 0.55
    )

    inf_prev <- estimate_inf_prev_from_model(pars$rho)

    stoch_80 <- sim_results() %>%
      filter(sigma==80) %>%
      select(rep, n_tested, n_recruited, n_fast, n_slow, n_overall)

    reps_i    <- unique(stoch_80$rep)
    cases_vec <- sapply(reps_i, function(r) estimate_case_incidence_from_model(stoch_80, pars=pars, my_rep=r))
    ARTI_vec  <- sapply(reps_i, function(r) estimate_ARTI_from_model(stoch_80, pars=pars, my_rep=r))

    tibble(
      Metric = c(
        "Infection prevalence, age 20",
        "Infection prevalence, age 30",
        "Infection prevalence, age 40",
        "ARTI (median)",
        "Cases per 100k per year (median)"
      ),
      Value = c(
        sprintf("%.2f", inf_prev["inf_prev_age20"]),
        sprintf("%.2f", inf_prev["inf_prev_age30"]),
        sprintf("%.2f", inf_prev["inf_prev_age40"]),
        sprintf("%.3f", median(ARTI_vec)),
        sprintf("%.0f", median(cases_vec))
      )
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

}


shinyApp(ui, server)
