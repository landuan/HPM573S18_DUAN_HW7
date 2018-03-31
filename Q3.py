import scipy.stats as Stat

##50% of the patients in our simulated model survived beyond 5 years:
percentage = 0.5

#a clinical study reports 400 of 573 participants survived:
kstudy=400
nstudy=573

#caculate the likelihood:
weight=Stat.binom.pmf(k=kstudy,n=nstudy,p=percentage)

print(weight)