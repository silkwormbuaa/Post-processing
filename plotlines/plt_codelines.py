#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   plt_codelines.py
@Time    :   2023/11/02 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   None
'''

import os
import sys
import matplotlib
import matplotlib.pyplot  as     plt
import matplotlib.dates   as     mdates
from   datetime           import datetime

plt.rcParams["text.usetex"] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size']   = 20

months = ["2022-12","2023-05","2023-08","2023-10","2023-11"]

lines = [6489, 13712, 18022,22558,24876]

date_objects = [datetime.strptime(month, "%Y-%m") for month in months]

fig, ax = plt.subplots(figsize=(10, 8))


ax.plot(date_objects, lines, marker='o', linestyle='--', color='black')

ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

plt.xticks(rotation=45)

fig.subplots_adjust(left=0.2, right=1, bottom=0.15, top=0.95)

plt.show()