#%%
import datetime
import requests
import json
import pandas as pd
import matplotlib.pyplot as plt
import requests_cache

# %%
df = pd.read_csv('https://covid.ourworldindata.org/data/owid-covid-data.csv')
# %%
df
# %%
df.keys()
# %%
df.dtypes
# %%
df = df[df.date == '2021-05-26']
# %%
df = df.set_index('iso_code')
# %%
inc = df.new_cases_smoothed_per_million / 10 * 7

# %%
df
# %%
df.new_cases_smoothed_per_million
# %%
