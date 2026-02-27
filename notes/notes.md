# 6 June 2025

Personal recap of the meeting this week with Kristin and Tyler:

-   **Aim 1:** How might a test for recency of TB infection (e.g., TASA) improve our ability to do the things listed below, beyond what's possible with existing diagnostics (e.g., IGRA)?
    -   Assessing risk of infection
    -   Measuring incidence/prevalence (more accurately, with less lag)
    -   Vaccine trials (measuring community-level transmission)
    -   Compare to Styblo approach
    -   Can we specify the sensitivity, specificity, and timing of a recency test to be maximally useful for these applications?
-   **Aim 2:** How might a test for TB infectiousness (CASS, facemasks) improve our ability to control TB?
    -   There's confounding between supershedders and super-contacters. How can we account for this? In which contact contexts might an indicator of biological infectiousness be helpful?
    -   Tools like CASS have been sidelined because they're not "effective" enough... but could this be a feature, not a bug? if they can only detect people with extremely high infectiousness, could this be exactly what we want?
    -   What kind of sensitivity and specificity -- in terms of predicting the number of secondary cases -- would we want from an infectiousness test? Can we specify the test parameters we'd want for it to be a useful outbreak control tool?

I'd like to code up some initial sims. I want to see if I've got the right ideas in mind, and we might be able to use some of the output for preliminary data.

Some code architecture:

-   Simulate an epidemic curve
-   Simulate sampling at various points in time, with different tests
-   Show estimates of incidence and prevalence over time

Before diving in, I want to look at existing TB models. Some useful resources:

-   [Guidance for country-level TB modelling](https://researchonline.lshtm.ac.uk/id/eprint/4653000/1/gomez_etal_2019_guidance_for_country-level_tb_modelling.pdf) (led by Nick Menzies)
-   [Progression from latent infection to active disease in dynamic tuberculosis transmission models: a systematic review of the validity of modelling assumptions](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(18)30134-8/abstract) (Menzies ... Cohen)
-   [Prospects for Tuberculosis Elimination in the United States: Results of a Transmission Dynamic Model](https://academic.oup.com/aje/article-abstract/187/9/2011/4995883) (Menzies, Cohen ... Salomon)
-   [Comparative Modeling of Tuberculosis Epidemiology and Policy Outcomes in California](https://www.atsjournals.org/doi/full/10.1164/rccm.201907-1289OC) (Menzies ... Shete)
-   [The Impact of Realistic Age Structure in Simple Models of Tuberculosis Transmission](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008479) (Brooks-Pollock, Cohen, Murray)
-   [Interferon-Gamma Release Assays versus Tuberculin Skin Testing for the Diagnosis of Latent Tuberculosis Infection: An Overview of the Evidence](https://onlinelibrary.wiley.com/doi/10.1155/2013/601737) (Trajman, Steffen, Menzies)

# 16 July 2025

Aim 3 -- if you could collect some bare minimum info on contacts, would that be enough to disentangle timing of infectiousness(?)

Some questions:

-   We now have two aims pages -- just tests for recent infection, or both tests for recency and tests for infectiousness?
    -   if we do scope it down to two aims: are those substantive enough to sustain two full aims?
    -   probably yes -- let's maybe just keep it with the two aims.
-   A potential [collaborator](https://wikitia.com/wiki/Claudia_Denkinger)

Scheduling -- deadline Oct 12 (Sun), internal deadline Oct 6-7. Sep 29 for admin documents.

9-ish weeks between now and then;

Aim for push over next 5 weeks; most of writing in August.

Significance & Innovation (2 weeks) Approach (2 weeks)

Check for draft from Kristin in early Aug...

Aim 2 feels like the more challenging of the two?

Plot two incidence curves -- ideally with similar IGRA profiles but yielding different TASA curves, just to show that there's information there.

(maybe something too to deminstrate the value of starting TBT for people with a recent infection vs. an IGRA-positive infection)

one more application in the context of vax trials? --\> gating entry into vaccine trial based on people who are IGRA positive. Could gating by TASA be better? (1) is transmission intense enough and (2) at an individual level, could we use TASA to be more predictive for good eligibility?

(still structure aims as vaccine trials + epidemic dynamics)

Look out for rough schedule to get us to submission. Review, give comments.

# 27 July 2025

Some goals for preliminary data modelling:

-   Simulate a clinical trial, or do some sample size calculations for a trial, under different rates of progression to TB. Consider using a test to (a) determine eligibility for a trial and (b) refine sample size estimates by getting better notion of incidence in a community.

Let's start there. A second aim would be similar to what I was working on above: how might a test for recency of infection reveal changes in incidence/prevalence that would be obscured with an IGRA-style test?

Let's try to come up with something simple for the first case:

For now, imagine TB incidence is at a steady state. Can we simulate how long ago a person was infected, given that we're at this steady state? And then we can simulate timing of infection, test positivity, etc.

Imagine $\rho$ is the incidence of TB (e.g. new cases per 100,000 people) in a given country. There's also going to be some rate of progressing to active TB after infection, and this I think is non-constant: you have a high risk close to your time of infection, then after that, the risk goes down. I should try to find some curves in the literature on this.

If I did have it in hand, what would I do? Maybe (a) simulate infection times, (b) simulate whether a person goes on to active disease (maybe that's the way to do it), (c) simulate the time the person goes on to active disease, for those who do; where I'd reduce the probability of progressing to active disease for those who are vaccinated.

(If I'm remembering correctly, I think that a main use case for vaccines is to prevent the progression to active disease; that's why we'd be recruiting people who are IGRA positive but aren't currently showing symptoms.)

Kristin notes the interesting fact that, using a test for recency of infection, it'd maybe take longer to recruit but less time to reach the endpoints. That said, since we'd need a smaller sample size, maybe we wouldn't actually need longer to recruit?

I'm going to work with Structure F from here: <https://pubmed.ncbi.nlm.nih.gov/29653698/>

Some useful information in the first paragraph of that paper:

> Latent infection is a defining feature of tuberculosis epidemiology. On infection with Mycobacterium tuberculosis, approximately 5% of otherwise healthy adults will develop active disease within 2 years (so-called fast progressors).1,2 Individuals who do not have rapid progression are classified as having slow-progressing latent tuberculosis infection. With latent infection, individuals experience no adverse health effects and will not transmit M tuberculosis, but they face an ongoing risk of developing active tuberculosis through reactivation. For individuals with long-established infection, the annual risk of active tuberculosis is low; empirical estimates are on the order of 10–20 per 100 000 individuals.3 However, as a result of high prevalence of latent tuberculosis infection in many settings,4reactivation can represent a substantial proportion of incident tuberculosis cases, or even the majority of such cases in settings in which transmission has been in sustained decline.5 The risk of progressing to active disease also varies by individual characteristics, with infants,6 individuals with advanced HIV infection,7,8 and individuals with other conditions that affect immune function9–12 having elevated progression risks.

So: aim for 5% of people as fast progressors who develop disease on average after 1 year. I might model this using an exponential distribution with rate 1.5 / year, which gives us 77% of people progressing after 1 year and 95% of people progressing after 2 years. We can play around with this, potentially increasing the rate if that's still not high enough, or making it a non-uniform hazard.

For the remaining 95% of people, we might give them a rate of .0001 of progressing to active disease.

For incidence, maybe we can use something like between 100 and 250 cases per 100,000 per year, representing a middle- to high-burden country (see page 9 here: <https://iris.who.int/bitstream/handle/10665/379339/9789240101531-eng.pdf?sequence=1>)

Then, we need to simulate test positivity: maybe we just do something deterministic, where both tests immediately turn positive (or do so after some short lag), and IGRA stays on always, while TASA drops off after, say 2 years. This is maybe also a parameter we'll want to vary.

# 30 July 2025

I did some useful work on the probability expressions to figure out the probability a person is infected and of a specific progression type, given that they're asymptomatic. That will be helpful, I think; but I also think that a brute-force algorithm could be useful, and ultimately will be more flexible. That's what I want to think through now.

The idea: imagine drawing a single random person from the population, one at a time. For each person, I'll draw:

-   their age
-   their progression type (fast or slow)
-   the timing of their infection (which could exceed their age, in which case the person isn't infected at the time of sampling)
-   the time at which they develop symptoms

I want to track:

-   Do we test the person? (We do this if they're asymptomatic; symptomatic people don't get tested and are rejected from the study outright)
-   Do we recruit the person into the study? (we do this if they're asymptomatic but test positive, i.e. if they're asymptomatic and infected within the last $\sigma$ years -- where $\sigma$ is equal to their age for IGRA, and might be something like 2 years for TASA)
-   How many people do we test? How many people do we recruit? What's the fraction of these?

Also: I need to think clearly through how we're doing the ultimate inference on vaccine ffectiveness. There are some issues here with a naive approach: if we follow up for five years and compare the number of progressions in the vaccinated vs. unvaccinated group, we'll get no difference, because even if vaccination makes it so that you take twice as long to progress, the fast progressors will all still generally progress within five years. Similarly -- if we just use a simple gamma exponential rate estimation with censoring, ther eare so many people who are censored that it completely washes out the information from the people who did progress. So, we might want something like a mixtur emodel, where we estimate whether a person is a slow or fast progressor, and then separately estimate their progression rate. This might make sense, because anyone who hasn't converted asfter five years is likely a slow progressor anyway, and they're not going to contribute much information -- and in fact we might not even care about them much, since slow progressors generally aren't as influential for disease transmission.

I wonder if a better mechanism of vaccine action might be to simply move people from the fast-progressor to the slow-progressor group with some probability equal to the VE.

This sort of works -- but really incidence needs to be higher if we're going to get the number of tests down to something like 8000

Good stuff in ncalc.R. Bed.

# 2 Aug 2025

Ok, so yesterday I got the derivation finally working in which I estimated the probability that a person is (a) infected in the past $\sigma$ years AND (b) a fast (or slow) progressor, given that they're currently symptomatic at age $a$.

This lets us simulate trial recruitment straightforwardly: if we test an asymptomatic person, this gives us the probability that they test positive on a test that detects infections within the past $\sigma$ years, and it tells us if the person is a fast or slow progressor. Of course, in reality, we won't know who is a fast or slow progressor, but knowing that for the simulations is important so that we cna project when the recruited people develop symptoms.

Now, I want to do something a little more thoughtful: given an age distribution, can we simulate recruitment? Here's the idea:

-   Draw a person of age $a$ from the population's age distribution, possibly restricting to [18, 50) to align with other trials
-   Given that person's age, calculate the joint probability that they're (a) infected in the past $\sigma$ years, (b) asymptomatic, and (c) a fast (slow) progressor. I think we should be able to do this using quantities I've already derived.

The reason we want this last thing is because, then, given a person's age, we can estimate the probability that we (a) test them and (b) what the outcome of that test is.

No, better: what we should do is

1.  Draw a person from the population, using the population's age distribution
2.  Determine if that person is asymptomatic, for which we can use the denominator of the expression I derived yesterday
3.  If they're asymptomatic, then we can assume we test them. Then, we can use the full expression I derived yesterday to calculate the probability that they test positive (infected in the last $\sigma$ years) and are a fast/slow progressor, given that they're asymptomatic.

I think that's the way forward.

Welp, after all that work... the brute force algorithm runs faster. No idea why. Maybe the thing to do is to stick with that and use the theory for intuition/plotting mean lines.

# 3 Aug 2025

Would love to get some decent simulations done today.

Excellent -- got the simulations done, and a theoretical curve plotted over the top that matches well. I think this is good enough to share.

# 19 Feb 2026 (LH)

Added basic README and example workflow. New naming conventions introduced: analytic and stoch (formerly theory and trials). Parameter data taken from specific historic trials is now also labelled by the vaccine prototype and phase number e.g. m72IIb (formerly nejm). Main running code is now run_analysis.qmd.

# 23/24 Feb 2026 (LH)

Age structure added to model - population age distribution and incidence by age. Works for "uniform" or any UN-recognised country.

# 25/26 Feb 2026 (LH)

Dealing with varying hazard rates for time to infection (sim_stoch() function). Initially we use an exponential distribution with constant annual incidence rho to estimate time-to-infection. Now rho is age-varying, so I've been investigating survival analysis and other ways to estimate 'time to event X'.

In short, my conclusion is to use a general survival function with empirically-specified hazard rate (inc_by_age). The general survival function collapses to exponential when the hazard rate is constant (which is what we want).

-   S(a) = exp(−∫₀ᵃ ρ(u) du) for age-specific incidence rho

-   This expression collapses to the exponential survival function when rho is constant

-   Discrete version: ∫₀ᵃ ρ(u) du = Σ (ρ(aᵢ) · Δaᵢ) i.e. adding up the cumulative incidence up to age a, for discrete age bands of width Δaᵢ.

See <https://en.wikipedia.org/wiki/Failure_rate#Conversion_to_cumulative_failure_rate> for a derivation of the equation relating S and rho (comes from solving a differential eqn).

From this survival function, we then compute the survival probability mass function (successive differences) and append a probability for tinf \> 100 (i.e. never infected during lifetime). One then draws a random number tinf from this sample, and procedure continues as before.

How has this changed model outputs? The distribution of tinf times now follows the survival pmf curve (and generally the inc_by_age curve):

![old_tinfs_recruited_histogram](images/old_tinfs_recruited_histogram-01.png)

![new_tinfs_recruited_histogram.png](images/new_tinfs_recruited_histogram.png)

Important note: Currently, incidence of *cases* is used as the hazard rate for tinf. What we actually want is incidence of *infection.*

# 27 Feb 2026 (LH)

Time to infection (tinf in sim_stoch()) currently uses incidence of *cases* as the hazard rate. We could instead use this hazard rate to calculate tsymp, and then backcalculate tinf values. This is tricky for tsymp \> 100 - I instead tried drawing a rexp(1, 0.0025) if the tsymp \> 100. This is a little better, gives a longer tail, but we still end up with negative tinf values. See plots here of 1000 draws of progressor_type and tsymp, with tinf then backcalculated:

![](images/calculating_tsymp_first---tsymps_drawn_histogram.png)

![](images/calculating_tsymp_first---tinfs_drawn_histogram.png)

Most tinfs are negative with this method.... Need to find another workaround. I will undo these changes.
