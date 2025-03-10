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
    data: input signal\n
    fs: sampling frequency\n
    n_seg: number of segments\n
    overlap: overlap ratio\n
    nfft: int, Length of the FFT used, if a zero padded FFT is desired. 
    If None, the FFT length is nperseg. Defaults to None.\n
    
    return: freqency array, psd
    """        
    
    len_perseg = int( len(data)//(1+(n_seg-1)*(1-overlap)) )
    
    len_overlap = int( len_perseg*overlap) 
    
    freq, psd = signal.welch( data, fs=fs, window='hann', 
                              nperseg=len_perseg, noverlap=len_overlap, 
                              nfft=nfft, scaling='density' )       
    
    return freq, psd

def pre_multi_psd( data, fs, n_seg, overlap, nfft=None, normalize=True ):
    
    """
    data: input signal\n
    fs: sampling frequency\n
    n_seg: number of segments\n
    overlap: overlap ratio\n
    nfft: int, Length of the FFT used, if a zero padded FFT is desired. 
    If None, the FFT length is nperseg. Defaults to None.\n
    
    return: freqency array, psd
    """        
    
    freq, psd = psd_hann( data, fs, n_seg, overlap, nfft=nfft )
    
    if normalize:
        power = np.sum( psd * freq[1] ) # frequency[1] is the frequency resolution
    else:
        power = 1.0
    
    
    pm_psd = freq * psd / power
    
    return freq, pm_psd

