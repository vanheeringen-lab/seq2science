import time

import seq2science


def log_welcome(logger, workflow):
    ascii_logo = (
f"""\
             ____  ____   __              
            / ___)(  __) /  \             
            \___ \ ) _) (  O )            
            (____/(____) \__\)            
                   ____                   
                  (___ \                  
                   / __/                  
                  (____)                  
   ____   ___  __  ____  __ _   ___  ____ 
  / ___) / __)(  )(  __)(  ( \ / __)(  __)
  \___ \( (__  )(  ) _) /    /( (__  ) _) 
  (____/ \___)(__)(____)\_)__) \___)(____)

version: {seq2science.__version__}
docs: https://vanheeringen-lab.github.io/seq2science
workflow: {workflow}
"""
    )

    logger.info(ascii_logo)

    # give people a second to appreciate this beautiful ascii art
    time.sleep(1)
