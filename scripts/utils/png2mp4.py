#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   png2mp4.py
@Time    :   2024/02/21 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   script of converting png files to mp4
'''

import os
import sys
import cv2

#https://stackoverflow.com/questions/63829991/qt-qpa-plugin-could-not-load-the-qt-platform-plugin-xcb-in-even-though-it
os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist

figpath = '/home/wencanwu/my_simulation/temp/220927_lowRe/snapshots/video_test/snapshots/figures_DS'

os.chdir( figpath )

figfiles = get_filelist( figpath, 'png' )

img = []

for i in range(len(figfiles)):
    img.append(cv2.imread(figfiles[i]))

height, width, layers = img[1].shape

# choose codec according to format needed
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
video = cv2.VideoWriter('video.mp4', fourcc, 5, (width, height))

for j in range(len(img)):
    video.write(img[j])

cv2.destroyAllWindows()
video.release()