# -*- coding: utf-8 -*-
'''
@File    :   vista_PSD.py
@Time    :   2022/11/10 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   For plotting power spectrum density
'''

import numpy             as     np
from   scipy             import signal


def psd_hann( data, fs, n_seg, overlap, nfft=None ):
    
    """
    data: input signal
    fs: sampling frequency
    n_seg: number of segments
    overlap: overlap ratio
    nfft: int, Length of the FFT used, if a zero padded FFT is desired. 
          If None, the FFT length is nperseg. Defaults to None.
    
    return: freqency array, psd
    """        
    
    len_perseg = int( len(data)//(1+(n_seg-1)*(1-overlap)) )
    
    len_overlap = int( len_perseg*overlap) 
    
    freq, psd = signal.welch( data, fs=fs, window='hann', 
                              nperseg=len_perseg, noverlap=len_overlap, 
                              nfft=nfft, scaling='density' )       
    
    return freq, psd

def pre_multi_psd( data, fs, n_seg, overlap, nfft=None ):
    
    """
    data: input signal
    fs: sampling frequency
    n_seg: number of segments
    overlap: overlap ratio
    nfft: int, Length of the FFT used, if a zero padded FFT is desired. 
          If None, the FFT length is nperseg. Defaults to None.
    
    return: freqency array, psd
    """        
    
    freq, psd = psd_hann( data, fs, n_seg, overlap, nfft=nfft )
    
    power = np.sum( psd * freq[1] ) # frequency[1] is the frequency resolution
    
    pm_psd = freq * psd / power
    
    return freq, pm_psd

