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
import gc
import sys
import cv2
from   moviepy.editor           import VideoFileClip, concatenate_videoclips

#https://stackoverflow.com/questions/63829991/qt-qpa-plugin-could-not-load-the-qt-platform-plugin-xcb-in-even-though-it
os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.tools       import get_filelist


def main():
    
    figpath = '/media/wencan/Expansion/temp/231124/postprocess/cf_wall/figs'
#    figpath = '/home/wencanwu/my_simulation/temp/220927_lowRe/snapshots/video_test/snapshots/figures_DS'

    os.chdir( figpath )

    figfiles = get_filelist( figpath, 'png' )
    
    n_fig = len(figfiles)
    size_figs = n_fig*os.path.getsize(figfiles[0])/1e9
    size_limit = 0.15
    
    if size_figs > size_limit:
        
        print(f"Total size of the images is {size_figs:10.3f} GB, \
                so we have to get video separately.")
        
        # number of videos and number of frames per video
        n_videos = int(size_figs//size_limit) + 1
        nfpv     = n_fig//n_videos
        
        print(f"The video is separated into {n_videos} parts.\n")
        
        # set the video names and the corresponding figfiles
        videos = [ f'output_{i:03d}.mp4' for i in range(n_videos) ]
        figs = [ figfiles[i*nfpv:(i+1)*nfpv] for i in range(n_videos-1) ]
        figs.append( figfiles[(n_videos-1)*nfpv:] )
        
        for i, (figfile, video) in enumerate(zip(figs, videos)):
            create_video_from_images(figfile, video, fps=25)
            print(f'total loading progress:{(i+1)/len(figs)*100:8.2f}%...')
        
        concatenate_videos(videos, '../output.mp4')

    else:
        create_video_from_images(figfiles, '../output.mp4', fps=25)
        
# ----------------------------------------------------------------------
# >>> Function Name                                                (Nr.)
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

def create_video_from_images(image_paths, output_video_path, fps=25):
    
    imgs = []
    
    for i, image in enumerate(image_paths):
        
        img = cv2.imread(image)
        
        # resize the image if it is too large
        height, width, layers = img.shape
        if width > 1920:
            img = cv2.resize(img, (1920, 1080))
            
        imgs.append(img)
        if i%10 == 0: print(f"{i} image loaded.")
    
    # Ensure all images have the same dimensions
    height, width, layers = imgs[0].shape
    
    print(f"Image dimensions: {height}x{width}x{layers}")
    
    for i in range(1, len(imgs)):
        if imgs[i].shape != (height, width, layers):
            print(f"Error: Image {image_paths[i]} has different dimensions.")
            return
    
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video  = cv2.VideoWriter(output_video_path, fourcc, fps, (width,height)) 
    
    for i in range(len(imgs)):
        video.write(imgs[i])
        if i%(len(imgs)//10) == 0: print(f"Progress : {i/len(imgs)*100:8.2f}%... ")
    
    cv2.destroyAllWindows()
    video.release()
    
    del imgs; gc.collect()


# ----------------------------------------------------------------------
# >>> Concatenate videos
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

def concatenate_videos(video_files, output_video_path):
    
    video_clips = []
    
    for video_file in video_files:
        video_clips.append(VideoFileClip(video_file))
    
    final_clip = concatenate_videoclips(video_clips)
    final_clip.write_videofile(output_video_path)
    
    for video_file in video_files:
        os.remove(video_file)


# ----------------------------------------------------------------------
# >>> Main()
# ----------------------------------------------------------------------
#
# Wencan Wu : w.wu-3@tudelft.nl
#
# History
#
# 2024/03/01  - created
#
# Desc
#
# ----------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())