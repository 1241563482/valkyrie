import pkg_resources, sys
sys.path.append('../')
from set_up import *

__version__ = "0.1"
__picture__ = f"""
 Version {__version__.ljust(5)}                                    
                 _  _                  _          
  /\   /\  __ _ | || | __ _   _  _ __ (_)  ___    
  \ \ / / / _` || || |/ /| | | || '__|| | / _ \   
   \ V / | (_| || ||   < | |_| || |   | ||  __/   
    \_/   \__,_||_||_|\_\ \__, ||_|   |_| \___|   
                           |___/                  
                                     by Yijie Zhu 
"""
__shell__ = pkg_resources.resource_filename('valkyrie', 'shell_scripts')
__python__ = pkg_resources.resource_filename('valkyrie', 'python_scripts')
