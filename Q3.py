import scipy.stats as Stat

# 50% of the patients in the model survived 5yr +
Percent = 0.5

# a clinical study reports 400 of 573 participants survived:
survived =400
total_studied=573

#caculate the likelihood:
weight=Stat.binom.pmf(k=survived,n=total_studied,p=Percent)

print(weight)