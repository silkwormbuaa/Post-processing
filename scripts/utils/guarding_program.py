#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@File    :   guarding_program.py
@Time    :   2024/11/07 
@Author  :   Wencan WU 
@Version :   1.0
@Email   :   w.wu-3@tudelft.nl
@Desc    :   check run.out file regularly and send email to user if the program is not running
          !!! this script should be soft linked to the case folder.
'''

import os
import sys
import time
from   sendgrid              import SendGridAPIClient
from   sendgrid.helpers.mail import Mail

source_dir = os.path.realpath(__file__).split('scripts')[0]
sys.path.append( source_dir )

from   vista.txt             import txt_tail

# =============================================================================

target_mail = 'w.wu-3@tudelft.nl'

# =============================================================================

# -----------------------------------------------------------------------------
# >>> read api key
def read_api_key():
    
    parameters = {}
    filename = os.path.join(os.path.expanduser("~"),".env")
    if not os.path.exists( filename ):
        raise FileExistsError(f"Please api key in {filename} first!")

    with open( filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.strip() and not line.startswith("#"):
                line = line.split("#")[0]
                key, value = line.strip().split("=",1)
                key = key.strip()
                value = value.strip().rstrip(",")
                parameters[key] = value
                
    return parameters

# -----------------------------------------------------------------------------

# =============================================================================
# >>> main program

# - read api key
keys = read_api_key()

api_key = keys['EMAIL_API_KEY']

# - get current working directory

casedir = os.path.basename( os.getcwd() )

casedir = "/media/wencanwu/Seagate Expansion Drive1/temp/test"

runfile = os.path.join( casedir, 'run.out' )

keep_guarding = True

time.sleep( 1200.0 )

while keep_guarding:

    if not os.path.exists(runfile):

        # when run.out file is not found, send an email to the user
        
        subject = f"Guarding Program: Error in {casedir}"
        plain_text_content = f"run.out file is not found in {casedir}!"

    else:
        
        # - check the timestamp of the run file.
        
        file_time    = os.path.getmtime( runfile )
        current_time = time.time()
        
        time_diff    = current_time - file_time
        
        # - when run.out file is not updated for more than 20 minutes, 
        #   check the content of run.out
        
        if time_diff > 1200.0:
        
            # - get the last 10 lines of the run file
            info = txt_tail( runfile )
        
            if "C'est tout! Danke und auf Wiedersehen." in info:
                
                subject = f"Guarding Program: Task landed!"
                plain_text_content = f"Task is finished safely in {casedir}!\n"
            
            elif "Expected time of arrival" in info:
                
                subject = f"Guarding Program: Task stopped strangely!"
                plain_text_content = f"Task stopped strangely in {casedir}! \n"
                
                
            keep_guarding = False
        
        else:
            
            # - if run.out file is just updated within 20 minutes, keep checking
            
            time.sleep( 1200.0 )

            
message = Mail( from_email='wuwencan@outlook.com',  # verified sender
                to_emails =target_mail,
                subject   =subject,
                plain_text_content=plain_text_content )

try:
    sg = SendGridAPIClient(api_key)
    response = sg.send(message)
    print(f"Email sent successfully with status code: {response.status_code}")
except Exception as e:
    print(f"Error: {e}")
