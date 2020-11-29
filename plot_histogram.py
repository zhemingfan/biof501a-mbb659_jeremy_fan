#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import argparse

df = pd.read_csv('results/lineage_report.csv', sep=',',index_col =0)

df['lineage'].value_counts().plot(kind='bar')

plt.ylabel('Frequency')
plt.xlabel('Lineage Types')
plt.title('Count of Different Lineage types within Samples')
plt.savefig('results/covid_histogram.png')