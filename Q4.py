from enum import Enum
import numpy as np
import scipy.stats as Stat
import scr.InOutFunctions as InOutSupport
import scr.StatisticalClasses as StatSupport
import scr.FormatFunctions as FormatSupport

N_TIME_STEPS=1000
ALPHA=0.05
NUM_COHORTS=1000
POST_L, POST_U, POST_N = 0.05, 0.25, 1000
POP_SIZE=573

class HealthStat(Enum):
    """ health status of patients  """
    ALIVE = 1
    DEAD = 0


class Patient(object):
    def __init__(self, id, mortality_prob):
        """ initiates a patient
        :param id: ID of the patient
        :param mortality_prob: probability of death during a time-step (must be in [0,1])
        """
        self._id = id
        self._rnd = np.random       # random number generator for this patient
        self._rnd.seed(self._id)    # specifying the seed of random number generator for this patient

        self._mortalityProb = mortality_prob
        self._healthState = HealthStat.ALIVE  # assuming all patients are alive at the beginning
        self._survivalTime = 0

    def simulate(self, n_time_steps):
        """ simulate the patient over the specified simulation length """

        t = 0  # simulation current time

        # while the patient is alive and simulation length is not yet reached
        while self._healthState == HealthStat.ALIVE and t < n_time_steps:
            # determine if the patient will die during this time-step
            if self._rnd.sample() < self._mortalityProb:
                self._healthState = HealthStat.DEAD
                self._survivalTime = t + 1  # assuming deaths occurs at the end of this period

            # increment time
            t += 1

    def get_survival_time(self):
        """ returns the patient survival time """
        # return survival time only if the patient has died
        if self._healthState == HealthStat.DEAD:
            return self._survivalTime
        else:
            return None


class Cohort:
    def __init__(self, id, pop_size, mortality_prob):
        """ create a cohort of patients
        :param id: cohort ID
        :param pop_size: population size of this cohort
        :param mortality_prob: probability of death for each patient in this cohort over a time-step (must be in [0,1])
        """
        self._patients = []      # list of patients
        self._survivalTimes = []  # list to store survival time of each patient
        self._five=[]

        # populate the cohort
        for i in range(pop_size):
            # create a new patient (use id * pop_size + n as patient id)
            patient = Patient(id * pop_size + i, mortality_prob)
            # add the patient to the cohort
            self._patients.append(patient)

    def simulate(self, n_time_steps):
        """ simulate the cohort of patients over the specified number of time-steps
        :param n_time_steps: number of time steps to simulate the cohort
        """
        # simulate all patients
        for patient in self._patients:
            # simulate
            patient.simulate(n_time_steps)
            # record survival time
            value = patient.get_survival_time()
            if not (value is None):
                self._survivalTimes.append(value)
                if value>5:
                    self._five.append(1)
                else:
                    self._five.append(0)

    def get_ave_survival_time(self):
        """ returns the average survival time of patients in this cohort """
        return sum(self._survivalTimes)/len(self._survivalTimes)

    def get_percentage(self):
        return float(sum(self._five))/float(len(self._five))

class MultiCohort:
    """ simulates multiple cohorts with different parameters """

    def __init__(self, ids, pop_sizes, mortality_probs):
        """
        :param ids: a list of ids for cohorts to simulate
        :param pop_sizes: a list of population sizes of cohorts to simulate
        :param mortality_probs: a list of the mortality probabilities
        """
        self._ids = ids
        self._popSizes = pop_sizes
        self._mortalityProbs = mortality_probs

        self._survivalTimes = []      # two dimensional list of patient survival time from each simulated cohort
        self._meanSurvivalTimes = []   # list of mean patient survival time for each simulated cohort
        self.percentages=[]

    def simulate(self, n_time_steps):
        """ simulates all cohorts """

        for i in range(len(self._ids)):
            # create a cohort
            cohort = Cohort(self._ids[i], self._popSizes[i], self._mortalityProbs[i])
            # simulate the cohort
            cohort.simulate(n_time_steps)
            # store average survival time for this cohort
            self._meanSurvivalTimes.append(cohort.get_ave_survival_time())
            # store percentages for this cohort
            self.percentages.append(cohort.get_percentage())

    def get_cohort_mean_survival(self, cohort_index):
        """ returns the mean survival time of an specified cohort
        :param cohort_index: integer over [0, 1, ...] corresponding to the 1st, 2ndm ... simulated cohort
        """
        return self._meanSurvivalTimes[cohort_index]


class Calibration:
    def __init__(self):
        self.cohortids=range(POST_N)
        self.mortalitysamples=[]
        self.mortalityresamples=[]
        self.weights=[]
        self.nomweights=[]
        self._csvRows = \
            [['Cohort ID', 'Likelihood Weights', 'Mortality Prob']]

    def sample_posterior(self):
        self.mortalitysamples = np.random.uniform(low=POST_L, high=POST_U, size=POST_N)

        multiCohort = MultiCohort(
            ids=self.cohortids,
            mortality_probs=self.mortalitysamples,
            pop_sizes=[POP_SIZE]*POST_N
        )
        ##simulate multi-cohort
        multiCohort.simulate(N_TIME_STEPS)
        ###caculate the likelihood
        self.weights=Stat.binom.pmf(k=400,n=573,p=multiCohort.percentages)
        ##caculate the weights
        self.nomweights=np.divide(self.weights,np.sum(self.weights))
        # re-sample mortality probability (with replacement) according to likelihood weights
        self.mortalityresamples = np.random.choice(
            a=self.mortalitysamples,
            size=NUM_COHORTS,
            replace=True,
            p=self.nomweights)

        # produce the list to report the results
        for i in range(0, len(self.mortalitysamples)):
            self._csvRows.append(
                [self.cohortids[i], self.nomweights[i], self.mortalitysamples[i]])

    def get_mortality_estimate_credible_interval(self, alpha, deci):
        """
        :param alpha: the significance level
        :param deci: decimal places
        :returns text in the form of 'mean (lower, upper)' of the posterior distribution"""

        # calculate the credible interval
        sum_stat = StatSupport.SummaryStat('Posterior samples', self.mortalityresamples)

        estimate = sum_stat.get_mean()  # estimated mortality probability
        credible_interval = sum_stat.get_PI(alpha)  # credible interval

        return FormatSupport.format_estimate_interval(estimate, credible_interval, deci)


m=Calibration()
m.sample_posterior()
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1-0.05, prec=0),
      m.get_mortality_estimate_credible_interval(0.05, 4))