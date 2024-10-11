# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 19:04:17 2018
    This class is a context manager which shows the running time of a code block.
@author: weibo
"""

import sys
import time

class timer():

    def __init__(self, message):
        self.message = message
        self.start = time.time()

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        # self.end = time.time()
        self.interval = time.time() - self.start
        hours, remainder = divmod(self.interval, 3600)
        minutes, seconds = divmod(remainder, 60)
        print("%s took: %02d:%02d:%02d" % (self.message, hours, minutes, seconds))
        print("=====================\n")

    def remainder(self, current_progress):
        
        elapse = time.time() - self.start
        lefttime = elapse * (1 - current_progress) / current_progress
        hours, _         = divmod(lefttime, 3600)
        minutes, seconds = divmod(_, 60)
        
        output = "Time left: %02d:%02d:%02d."%(hours, minutes, seconds)
        
        return output
    
    def print_progress( self, i, size, rank=None):
        progress = (i+1)/size
        print(f"Rank {rank}: {i+1}/{size} is done. " + self.remainder(progress))
        print("------------------\n")
        sys.stdout.flush()
