sim_stoch <- function(pars, fasttarget=50, agedist){
  with(as.list(pars), {
    
    # Initialize tracking variables  
    n_tested <- 0
    n_recruited <- 0
    n_fast <- 0
    n_overall <- 0
    recruited_df <- tibble()
    
    while(n_fast < fasttarget){
      # Grab their age from the age distribution
      age <- sample(minage:maxage, size=1, prob=agedist[minage:maxage])  # R normalises the subset automatically
      
      # Grab their progression status 
      progressor_type <- sample(c("slow","fast"), 
                                size=1, prob=c(p_slow, 1-p_slow))
      
      # Simulate their time of infection
      tinf <- rexp(1, rho)
      
      # Simulate their time of symptoms 
      tsymp <- tinf + rexp(1, (if(progressor_type=="slow"){mu_slow} else {mu_fast}))
      
      # ELIGIBILITY
      if(tsymp > age){
        # They are asymptomatic, so let's test them: 
        n_tested <- n_tested + 1
        if((tinf >= age-sigma) & (tinf <= age)){
          # They were infected in past sigma years, so test positive. Recruit!
          n_recruited <- n_recruited + 1
          if(progressor_type=="fast"){n_fast <- n_fast + 1}
          recruited_df <- bind_rows(recruited_df, tibble(  # bind_rows inside a loop is slow
            id=n_recruited,
            age=age,
            progressor_type=progressor_type,
            tinf=tinf
          ))
        }
      }
      
      n_overall <- n_overall + 1
    }
    
    out <- list(n_tested=n_tested, recruited_df=recruited_df)
    
    return(out)
  })
}  # 