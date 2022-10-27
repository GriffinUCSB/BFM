import pandas as pd
import numpy as np
import sys

input = "merged_all_genes.csv"

df = pd.read_csv(input)

rank = []
for i in range (len(df)):
    rank.append(i)

#rank ds epsilon values
df = df.sort_values(by="ds epsilon")
df["rank dse"] = rank

#rank ss epsilon values
df = df.sort_values(by="ss epsilon")
df["rank sse"] = rank

#compare ranks
df["rank-order comparison"] = df["rank dse"]-df["rank sse"]

#median transform ds epsilon and ss epsilon so std ~ 1 and mean ~ 0
ds_median = df["ds epsilon"].median()
df["ds median transformed epsilon"] = (df["ds epsilon"] - ds_median)/(df["ds epsilon"].std())
print("ds:: std = " + str(df["ds median transformed epsilon"].std()) + " ;mean = " + str(df["ds median transformed epsilon"].mean()) + " ;median = " + str(df["ds median transformed epsilon"].median()))

ss_median = df["ss epsilon"].median()
df["ss median transformed epsilon"] = (df["ss epsilon"] - ss_median)/(df["ss epsilon"].std())
print("ss:: std = " + str(df["ss median transformed epsilon"].std()) + " ;mean = " + str(df["ss median transformed epsilon"].mean()) + " ;median = " + str(df["ss median transformed epsilon"].median()))

df["median normalized comparisson"] = df["ds median transformed epsilon"]- df["ss median transformed epsilon"]
df.to_csv("rank_+_normalized_median_comparison.csv", sep=",", header=True, index=False, mode="a")
