# 6 June 2025 

Personal recap of the meeting this week with Kristin and Tyler: 

- **Aim 1:** How might a test for recency of TB infection (e.g., TASA) improve our ability to do the things listed below, beyond what's possible with existing diagnostics (e.g., IGRA)? 
	- Assessing risk of infection
	- Measuring incidence/prevalence (more accurately, with less lag) 
	- Vaccine trials (measuring community-level transmission) 
	- Compare to Styblo approach 
	- Can we specify the sensitivity, specificity, and timing of a recency test to be maximally useful for these applications? 
- **Aim 2:** How might a test for TB infectiousness (CASS, facemasks) improve our ability to control TB? 
	- There's confounding between supershedders and super-contacters. How can we account for this? In which contact contexts might an indicator of biological infectiousness be helpful? 
	- Tools like CASS have been sidelined because they're not "effective" enough... but could this be a feature, not a bug? if they can only detect people with extremely high infectiousness, could this be exactly what we want? 
	- What kind of sensitivity and specificity -- in terms of predicting the number of secondary cases -- would we want from an infectiousness test? Can we specify the test parameters we'd want for it to be a useful outbreak control tool? 


I'd like to code up some initial sims. I want to see if I've got the right ideas in mind, and we might be able to use some of the output for preliminary data. 

Some code architecture: 

- Simulate an epidemic curve 
- Simulate sampling at various points in time, with different tests 
- Show estimates of incidence and prevalence over time 

Before diving in, I want to look at existing TB models. Some useful resources: 

- [Guidance for country-level TB modelling](https://researchonline.lshtm.ac.uk/id/eprint/4653000/1/gomez_etal_2019_guidance_for_country-level_tb_modelling.pdf) (led by Nick Menzies) 
- [Progression from latent infection to active disease in dynamic tuberculosis transmission models: a systematic review of the validity of modelling assumptions](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(18)30134-8/abstract) (Menzies ... Cohen) 
- [Prospects for Tuberculosis Elimination in the United States: Results of a Transmission Dynamic Model](https://academic.oup.com/aje/article-abstract/187/9/2011/4995883) (Menzies, Cohen ... Salomon)
- [Comparative Modeling of Tuberculosis Epidemiology and Policy Outcomes in California](https://www.atsjournals.org/doi/full/10.1164/rccm.201907-1289OC) (Menzies ... Shete) 