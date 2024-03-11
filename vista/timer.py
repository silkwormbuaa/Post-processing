# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 19:04:17 2018
    This class is a context manager which shows the running time of a code block.
@author: weibo
"""

import time


class timer():

    def __init__(self, message):
        self.message = message

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


