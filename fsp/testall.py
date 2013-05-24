#!/usr/bin/env python

""" 
    Runs all test jobs for the SPAC smooth particle code.
    Renames and copies images and generates latex source for the section
    of the manual in which the test plots are recorded.

    Pass 'report_only' as a first argument to skip the tests
    and just rebuild the report. Of course this will only work
    if you have previously run the tests at some point.

"""

import os
import run
import report
import sys

#RUNTESTS = True
#REPORT = True

#if len(sys.argv) > 1:
#    if sys.argv[1] == 'report_only':
#        RUNTESTS = False

#if RUNTESTS:
for i in range(1,11):
        print "Running test",i
        os.system("./run.py --test " + str(i))

#if REPORT:
#    print "Generating test report"
#    print run.TEST_ODIR
    #TNAME = run.TEST_NAME

#    for key,name in TDIR.items():
#        print run.odir + '/' + name
#        break

    # The report-writing is implemented in vasp/report.py?
    # or is it in manual_report.py and report.py has been abandoned?


