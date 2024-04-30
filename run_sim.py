"""
Define the HPVsim simulation
"""

# Standard imports
import pylab as pl
import numpy as np
import sciris as sc
import hpvsim as hpv

# %% Settings and filepaths
# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True


# %% Simulation creation functions
def make_st(product='hpv', screen_coverage=0.15, treat_coverage=0.7, start_year=2020):
    """ Make screening & treatment intervention """

    # Define who's eligible for screening
    age_range = [30, 50]
    len_age_range = (age_range[1]-age_range[0])/2
    # Targeting x% of women aged 30-50 implies that the following proportion need to be screened each year
    model_annual_screen_prob = 1 - (1 - screen_coverage)**(1/len_age_range)

    # Women are eligible for screening if it's been at least 5 years since their last screen
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
                                  (sim.t > (sim.people.date_screened + 5 / sim['dt']))

    # Make the routine screening interventions
    screening = hpv.routine_screening(
        prob=model_annual_screen_prob,
        eligibility=screen_eligible,
        start_year=start_year,
        product=product,
        age_range=age_range,
        label='screening'
    )

    # People who screen positive are assigned treatment
    screen_positive = lambda sim: sim.get_intervention('screening').outcomes['positive']
    assign_treatment = hpv.routine_triage(
        start_year=start_year,
        prob=1.0,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner'
    )

    ablation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=treat_coverage,
        annual_prob=False,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation'
    )

    excision_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['excision'].tolist() +
                                             sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(
        prob=treat_coverage,
        annual_prob=False,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )

    radiation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=treat_coverage/4,  # assume an additional dropoff in CaTx coverage
        annual_prob=False,
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    st_intvs = [screening, assign_treatment, ablation, excision, radiation]

    return st_intvs


def make_sim(location='nigeria', calib_pars=None, debug=0, interventions=None, analyzers=None, datafile=None, seed=1, end=2020):
    """
    Define parameters, analyzers, and interventions for the simulation
    """

    # Basic parameters
    pars = sc.objdict(
        n_agents=[20e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        end=end,
        genotypes=[16, 18, 'hi5', 'ohr'],
        location=location,
        ms_agent_ratio=100,
        verbose=0.0,
        rand_seed=seed,
    )

    # Sexual behavior parameters
    # Debut: derived by fitting to 2018 DHS
    # Women:
    #           Age:   15,   18,   20,   22,   25
    #   Prop_active: 19.2, 57.3, 73.9, 85.5, 91.6
    # Men:
    #           Age:  15,   18,   20,   22,   25
    #   Prop_active: 3.1, 14.5, 30.1, 51.9, 70.1
    # For fitting, see https://www.researchsquare.com/article/rs-3074559/v1
    pars.debut = dict(
        f=dict(dist='lognormal', par1=17.41, par2=2.75),
        m=dict(dist='lognormal', par1=17.91, par2=2.83),
    )

    # Participation in marital and casual relationships
    # Derived to fit 2018 DHS data
    # For fitting, see https://www.researchsquare.com/article/rs-3074559/v1
    pars.layer_probs = dict(
        m=np.array([
            # Share of people of each age who are married
            [0, 5, 10,      15,     20,     25,     30,     35,     40,     45,   50,   55,   60,   65,    70,   75],
            [0, 0,  0,  0.1596, 0.4466, 0.5845, 0.6139, 0.6202, 0.6139, 0.5726, 0.35, 0.21, 0.14, 0.07, 0.035, 0.007],
            [0, 0,  0,   0.228,  0.638,  0.835,  0.877,  0.886,  0.877,  0.818,  0.5,  0.3,  0.2,  0.1,  0.05, 0.01 ],
        ]),
        c=np.array([
            # Share of people of each age in casual partnerships
            [0, 5,  10,  15,  20,  25,  30,  35,    40,   45,   50,   55,   60,   65,   70,   75],
            [0,  0, 0.2, 0.6, 0.6, 0.6, 0.4, 0.4,  0.4,  0.1, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02],
            [0,  0, 0.2, 0.6, 0.8, 0.6, 0.6, 0.8, 0.95, 0.95, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
        ])
    )

    pars.m_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.2),
    )
    pars.f_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.2),
    )

    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    # Interventions
    sim = hpv.Sim(pars=pars, interventions=interventions, analyzers=analyzers, datafile=datafile)

    return sim


# %% Simulation running functions
def run_sim(
        location=None, analyzers=None, interventions=None, debug=0, seed=1, verbose=0.2,
        do_save=False, end=2020, calib_pars=None, meta=None):

    # Make sim
    sim = make_sim(
        location=location,
        debug=debug,
        interventions=interventions,
        analyzers=analyzers,
        calib_pars=calib_pars,
        end=end,
    )
    sim['rand_seed'] = seed
    sim.label = f'{location}--{seed}'

    # Store metadata
    sim.meta = sc.objdict()
    if meta is not None:
        sim.meta = meta  # Copy over meta info
    else:
        sim.meta = sc.objdict()
    sim.meta.location = location  # Store location in an easy-to-access place

    # Run
    sim['verbose'] = verbose
    sim.run()
    sim.shrink()

    if do_save:
        sim.save(f'results/{location}.sim')

    return sim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    # Make a list of what to run, comment out anything you don't want to run
    to_run = [
        # 'run_single',
        'run_scenario',
    ]

    location = 'nigeria'
    # st_intv = make_st(screen_coverage=0.15, treat_coverage=0.7)  # Make screening & treatment interventions

    # Run and plot a single simulation
    # Takes <1min to run
    if 'run_single' in to_run:
        sim = run_sim(interventions=None, end=2020, debug=debug)  # Run the simulation
        sim.plot()  # Plot the simulation

    # Example of how to run a scenario with different screening coverage rates
    # Takes ~2min to run
    if 'run_scenario' in to_run:
        baseline_st = make_st(screen_coverage=0.05, treat_coverage=0.3)  # Make screening & treatment interventions
        scenario_st = make_st(screen_coverage=0.50, treat_coverage=0.7)  # Imrpove screening & treatment coverage

        sim_baseline = run_sim(end=2060, interventions=baseline_st, debug=debug)
        sim_scenario = run_sim(end=2060, interventions=scenario_st, debug=debug)

        # Now plot cancers under the two alternative screen & treat scenarios
        pl.figure()
        res0 = sim_baseline.results
        res1 = sim_scenario.results
        what = 'cancers'
        pl.plot(res0['year'][50:], res0[what][50:], label='Baseline')
        pl.plot(res0['year'][50:], res1[what][50:], color='r', label='Improved screen & treat')
        pl.legend()
        pl.title(what)
        pl.show()

    # To run more complex scenarios, you may want to set them up in a separate file

    T.toc('Done')  # Print out how long the run took

