#!/Work/pipeline/software/Base/miniconda3/bin/python
# -*- coding: utf-8 -*-

import re
import os
import sys
import time
import logging
import argparse

from selenium import webdriver
from selenium.webdriver.chrome.options import Options


LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def html2png(url, out, paths="", times=5):

    jenkinsJobName=os.getenv("JOB_NAME")

    driver = webdriver.Edge("/Work/pipeline/software/Base/phantomjs/msedgedriver")
 #   try:
    #    driver = webdriver.Edge(executable_path="/Work/pipeline/software/Base/phantomjs/msedgedriver")
        #wait = WebDriverWait(driver, 10)
        #driver = webdriver.Firefox(executable_path=str(paths)) #geckodriver
  #  except:
        #driver = webdriver.PhantomJS(executable_path=str(paths)) #phantomjs
  #      driver = webdriver.Chrome(executable_path="chromedriver") #There is a problem using the Chrome browser.
    driver.maximize_window()
    driver.get(url) # Load page
    time.sleep(times) 
    driver.save_screenshot(out)
    driver.close()


def add_hlep_args(parser):
    parser.add_argument('--url', metavar='FILE', type=str, required=True,
        help="Input the path of the html file.")
    parser.add_argument('-p', '--paths', metavar='FLOAT', type=str,
        default="/Work/pipeline/software/Base/geckodriver/v0.29.1/bin/geckodriver",
        help="Input the browser path used,default=/Work/pipeline/software/Base/geckodriver/v0.29.1/bin/geckodriver.")
    parser.add_argument('-t', '--times', metavar='FLOAT', type=float, default=30,
        help="The time it takes for the browser to load the file,The larger the file, the more time it takes to set (s),default=30.")
    parser.add_argument('-o' ,'--out', metavar='FILE', type=str, default="txt.png",
        help="Set output file name,default=txt.png")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    html2png.py Convert files in html format to other formats.

attention:
    #export PATH=/nextomics/Software/meta/phantomjs/v2.1.1/bin:$PATH
    html2png.py --url "txt.html" -o txt.png

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    html2png(args.url, args.out, args.paths, args.times)


if __name__ == "__main__":

    main()
