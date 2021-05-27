# %%
import datetime
import requests
import json
import pandas as pd
import matplotlib.pyplot as plt
import requests_cache

# %%
requests_cache.install_cache()
# %%
# %%
params = {'pangolin_lineage': 'P.1',
          'location_id': 'DEU', 'mutations': 'S:E484K'}
# params = {'location_id': 'DEU'}
r = requests.get(
    'https://api.outbreak.info/genomics/prevalence-by-location', params)
r
# %%
raw_result = json.loads(r.content)
raw_result['results']  # list of results
# turn results into dataframe
# date, total count, lineage count, then can cumulate

# %%
# Method 1 Rolling Mean
df = pd.DataFrame(raw_result['results'])
df.date = pd.to_datetime(df.date)
df = df.set_index('date')
# %%
index = pd.date_range('20200101', datetime.datetime.today())
df2 = pd.DataFrame({'total_count': 0}, index=index)
df2.update(df)
df2.plot()
# %%
df2.total_count.rolling(14).sum().plot()
# %%
# (df.lineage_count / df.total_count).plot()
linlast = df.lineage_count.rolling(14).sum()
totlast = df.total_count.rolling(14).sum()
plt.plot(linlast/totlast)
# %%
# Method 2 Exponentially weighted: but need all sequences as well
linewm = df.lineage_count.ewm(halflife=7).mean()
totewm = df.total_count.ewm(halflife=7).mean()
plt.plot(linewm/totewm)
# %%
'''
Plan:
1. Specify Lineage
2. Specify Country
3. Specify Recency cutoff
4. Get Proportion
'''

# %%
r = requests.get(
    'https://api.outbreak.info/genomics/lineage-by-sub-admin-most-recent')
raw_result = json.loads(r.content)
countries = pd.DataFrame(raw_result['results'])
countries.set_index('name', inplace=True)
codes = pd.DataFrame({'full_name': countries.index, 'id': countries.id})
codes.set_index('id', inplace=True)

codes.drop('None', inplace=True)
# %%
# %%


def get_tot_sequences(id, lineage='', mut='', time=35):
    print(id)
    params = {'location_id': id}
    if lineage:
        params['pangolin_lineage'] = lineage
    if mut:
        params['mutations'] = mut
    r = requests.get(
        'https://api.outbreak.info/genomics/prevalence-by-location', params)
    raw_result = json.loads(r.content)
    if len(raw_result['results']) == 0:
        return 0
    else:
        df = pd.DataFrame(raw_result['results'])
        df.date = pd.to_datetime(df.date)
        df.set_index(df.date, inplace=True)
        df.drop('date', 1, inplace=True)
        index = pd.date_range('20200101', datetime.datetime.today())
        df2 = pd.DataFrame({'lineage_count': 0}, index=index)
        df2.update(df)
        if time == -1:
            result = df2.lineage_count.sum()
        else:
            result = df2.lineage_count.rolling(time).sum().iloc[-1]
        print(result)
        return result


# %%
# %%
codes['total_count'] = codes.apply(
    lambda row: get_tot_sequences(row.name), axis=1)

# %%
codes['B1351'] = codes.apply(
    lambda row: get_tot_sequences(row.name, 'B.1.351') if row.total_count > 0 else 0, axis=1)
# %%
codes = codes.assign(B1351_ratio=(1+codes.B1351)/(1+codes.total_count))
codes.B1351_ratio.sort_values(ascending=False).iloc[0:20]
# %%
codes['P1'] = codes.apply(
    lambda row: get_tot_sequences(row.name, 'P.1') if row.total_count > 0 else 0, axis=1)
# %%
codes = codes.assign(P1_ratio=(0+codes.P1)/(0+codes.total_count))
codes.P1_ratio.sort_values(ascending=False).iloc[0:20]
# %%
codes['India'] = codes.apply(
    lambda row: get_tot_sequences(row.name, 'B.1.617.2', time=28) if row.total_count > 0 else 0, axis=1)
# %%
codes = codes.assign(India_ratio=(0+codes.India)/(0+codes.total_count))
codes.India_ratio.sort_values(ascending=False).iloc[0:20]
# %%
codes = codes.assign(India_ratio_corr=(1+codes.India)/(1+codes.total_count))
codes.India_ratio_corr.sort_values(ascending=False).iloc[0:20]
# %%
x = range(20)
y = []
for i in x:
    y.append(get_tot_sequences('BLR', 'P.1', i) /
             get_tot_sequences('BLR', '', i))
plt.plot(x, y)

# %%
get_tot_sequences('BLR', 'P.1')
# %%
codes.total_count = codes.total_count.astype('int')
codes.to_csv('total_sequences.csv')
# %%
df = pd.read_csv('https://covid.ourworldindata.org/data/owid-covid-data.csv')
df = df[df.date == '2021-05-26']
# %%
df = df.set_index('iso_code')
# %%
inc = df.new_cases_smoothed_per_million / 10 * 7
# %%
vinc = pd.concat([codes, inc.rename('inc')], axis=1)
# %%
vinc = vinc.assign(combined=vinc.P1_ratio * vinc.inc)
# %%
vinc[vinc.total_count > -
     1].sort_values(by='combined', ascending=False).to_csv('vinc.csv')
# %%
