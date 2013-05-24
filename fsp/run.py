#!/usr/bin/env python
""" Run control script for the SPAC smooth particle program.
    Manages the creation of run directories, archiving of the
    source used, and calling of the post-processing script.
    Can also be used to run tests, and to continue
    previous runs.

    Andrew Charles 2013 All rights reserved.
    School of Applied Science, RMIT University, Melbourne.
   
    Todo
    ----
    - Use comands module to get output status of all commands we call 
    - Implement preprun and moverun in python (more
      portable).

   ./run.py -n RUN_NAME -f CONFIG_NAME

"""

def help():
    print """ 
    Run control script for the SPAC smooth particle program.
    Manages the creation of run directories, archiving of the
    source used, and calling of the post-processing script.
    Can also be used to run tests, and to continue
    previous runs.

    Command Line Options
    --------------------
    -f <config_File_name> This is the configuration file to use
    -n  <run name> This is the run name.
    -c  <run name>  Continue this run. Run name is a folder name.
    -s  <state_name> Name of file that holds initial state
    -a Use a new config file for the continued run (Only if -c defined).
    -t  --test <test case number> Run this test case. (Use the script
        testall.py to run all tests.  
    -b --base <directory> Base directory for this run.
    """

import shutil
import sys
import os
import string
import glob
from fsp_post import post, inflate_sph_ncdf, load_sphvars
import getopt
import os
sys.path.append("./fsp_post")

TEST = False
NORMAL = True 
INPLACE = False
CONTINUE = False
NEW_CONF = False
SET_CONFIG = False
CONFIG_FILE = ''
SET_STATE = False
STATE_FILE = ''
RUN_DIR = "."
RUN_NAME = "default"
test_ref = 0
BASE_DIR = os.environ.get("FSP_BASE")
RUN_BASE = os.environ.get("FSP_RUN")
TEST_DIR = BASE_DIR + "/test_cases"
TEST_ODIR = TEST_DIR + "/test_output" # output directory for tests
START_DIR = os.getcwd()

class FsphRun:
    """ Holds all the information required to run an instance of the
        fortran smooth particle code fsph.
        input_file -- the name of the run setup file containing such 
            parameters as the number of timesteps and particles.
        configuration -- this means the spatial configuration of the
            particles.

    """
    def __init__(self,input='none',config='none',name='none',
        description='description',
        long_description='long_description',
        plot=post.post_process):
        self.input_file = input
        self.config_file = config
        self.name = name
        self.directory = TEST_ODIR + '/' + self.name
	if not os.path.exists(self.directory):
		os.makedirs(self.directory)
        self.plot_command = plot
        self.description='description'
        self.long_description='long_description'
        # Other things you could do here include checking that
        # paths and files actually exist

    def run(self):
        """ Ok this is not all that clean, mostly calls shell scripts
            to do the heavy lifting of running the test.
        """
        print 'Running test in ',self.directory
        os.system("./preprun " +  self.directory)
        shutil.copyfile(self.input_file,self.directory + "/spinput.txt")
        if self.config_file != "none":
            shutil.copyfile(self.config_file,
                self.directory + "/sphstate.input")
        os.system("./moverun " +  self.directory)
        os.chdir(self.directory)
        os.system("./fsph")
        self.plot_command(self.name)

test = {}
test[0] = FsphRun(input=TEST_DIR+'/test_input_0',config="none",name='test_zero')

test[1] = FsphRun(input=TEST_DIR + "/test_input_1",
                   name="test_droplet",
                   description="100 particle droplet formation",
                   long_description="Droplet formation"
                   )

test[2] = FsphRun(input=TEST_DIR + "/test_input_2",
                   name="test_collision",
                   description="100 particle wall collision test",
                   long_description="This tests collisions with wall \
                    boundaries and  \
                    between particles, and also provides a test of the \
                    models energy conservation."
                    )

test[3] = FsphRun(input=TEST_DIR + "/test_input_3",
                   config=TEST_DIR + "/test_config_3",
                   name="test_grad",
                   description= "gradient term film width test")


test[4] = FsphRun(input=TEST_DIR + "/test_input_4",
                   name="test_film",
                   description="256 particle film formation",
                   long_description="This tests the formation of a  \
                        liquid-vapour film \
                        in a long box.")

test[5] = FsphRun(input=TEST_DIR + "/test_input_5",
                   name="test_settle",
                   description="256 particle settling liquid",
                   long_description="Settle")

test[6] = FsphRun(input=TEST_DIR + "/test_input_6",
                   name="test_bounds",
                   description="256 particle fill box")

test[7] = FsphRun(input=TEST_DIR + "/test_input_lf",
                   name="test_thermostat",
                   description="256 particle expansion with thermostat",
                   long_description="A fast expansion of a gas in a  \
                        periodic box \
                        strenuously tests the energy conservation of the \
                        scaling thermostat.")

test[8] = FsphRun(input=TEST_DIR + "/test_input_lf",
                   name="test_lf",
                   description="leapfrog integrator test")

test[9] = FsphRun(input=TEST_DIR + "/test_input_ie",
                   name="test_ie",
                   description="improved euler integrator test",
                   long_description="Fill")

test[10] = FsphRun(input=TEST_DIR + "/test_input_rk",
                    name="test_rk",
                    description="runge kutta integrator test")

test[11] = FsphRun(input=TEST_DIR + "/test_input_adke",
                    name="test_adke",
                    description="adaptive kernel droplet formation test",
                    long_description = "Tests the gradient term effect of  \
                        film interface width")

# CALLING THIS test12 for now...
test[11] = FsphRun(input=TEST_DIR + "/test_input_adke",
                    name="test_adke",
                    description="adaptive kernel droplet formation test",
                    long_description = "Tests the gradient term effect of  \
                        film interface width")

            
def process_args():
    global NORMAL,TEST,INPLACE,RUN_DIR,CONTINUE,test_ref,RUN_NAME,RUN_BASE
    global prev, NEW_CONF, SET_CONFIG, CONFIG_FILE, SET_STATE, STATE_FILE

    if(len(sys.argv) == 1):
        print 'No run name specified,running in place'
        INPLACE = True
        RUN_DIR = "."
        RUN_NAME = ""
    else:
        RUN_NAME = sys.argv[1]
        print "Treating arg1 as directory name"

    try:
        opts,args = getopt.getopt(sys.argv[1:],"n:t:hc:b:af:s:", \
            ["name=","test=","help","base=","continue=","base=","alltest","config=","state="])
    except getopt.GetoptError:
        print "Command line error"
        sys.exit(2)

    for opt,arg in opts:

        if opt in ("-b","--base"):
            NORMAL = False
            RUN_BASE = arg

        if opt in ("-f","--config"):
            SET_CONFIG = True
            CONFIG_FILE = arg
            print 'using config file',arg

        if opt in ("-t","--test"):
            print "Running test case " + arg
            TEST = True
            test_ref = int(arg)
            #print TEST_DESCRIPTION[test_ref]
            NORMAL = False

        if opt in ("-n","--name"):
            print "Setting name to " + arg
            RUN_NAME = arg
            NORMAL = False

        if opt in ("-s","--state"):
            print "Using state from file " + arg
            SET_STATE = True
            STATE_FILE = arg
        
        if opt in ("-c","--continue"):
            print "Continuing ",arg
            CONTINUE = True
            NORMAL = False
            prev = arg
            
        if opt in ("-a"):
            print "New config"
            NEW_CONF = True

        if opt in ("-h","--help"):
            help()
            sys.exit()

        #if opt in ("--alltest"):
        #    ALLTEST = True
        # see script testall.py
             
    RUN_DIR = RUN_BASE + "/" + RUN_NAME

    print "Run directory:",RUN_DIR


def run_test(i):
    test[i].run()


def main():

    if os.path.isfile('fsph'):
        print 'Found SPH executable'
    else:
        print 'SPH executable not found'

    if INPLACE:
        # Run and render in place. Usually only done while working
        # on the code.
        print "Running in place"
        os.system("./fsph ")
        post.post_process('default')

    if NORMAL:
        # This is a normal run 
        print 'Starting a normal run.'
        print 'Creating directory '+ RUN_DIR
        os.system("./preprun " + RUN_DIR)
        os.chdir(RUN_DIR)
        print 'Running smoothed particle program...'
        os.system("./fsph ")
        post.post_process(RUN_NAME)

    elif SET_CONFIG:
        print 'Creating directory '+ RUN_DIR
        print 'Starting a preconf run.'
        
        if SET_STATE: 
            fname = RUN_BASE + "/" + STATE_FILE + "/" + STATE_FILE + ".nc"
            print 'using',fname
            if not (os.path.exists(RUN_DIR)):
	        os.mkdir(RUN_DIR)
            #in_name = RUN_DIR + "/sphstate.input"
            in_name = "sphstate.input"
            # Use inflate_sph_netcdf to extract last step from netcdf file
            inflate_sph_ncdf.output_step(fname,-1,in_name)
            
        os.system("./preprun " + RUN_DIR + ' ' + CONFIG_FILE)
        os.chdir(RUN_DIR)
        print 'Running smoothed particle program...'
        os.system("./fsph ")
        post.post_process(RUN_NAME)
        os.chdir(START_DIR)


    elif TEST:
        current_test = test[test_ref]
        current_test.run()

    else:
        print 'Creating directory '+ RUN_DIR
        os.system("./preprun "+ RUN_DIR)
        if CONTINUE:
            fname = RUN_BASE + "/" + prev + "/" + prev + ".nc"
            in_name = RUN_DIR + "/sphstate.input"

            # Use inflate_sph_netcdf to extract last step from netcdf file
            inflate_sph_ncdf.output_step(fname,-1,in_name)

            if NEW_CONF:
                shutil.copyfile("spinput.txt",RUN_DIR+"/spinput.txt")            
            else:
                # Use the old conf file but modify the load state att
                config_name = RUN_BASE + "/" + prev + "/spinput.txt"
                # Modify spinput to use the state of the previous run.
                # Load the config into a dictionary
                config = load_sphvars.load_sphvars_new(open(config_name,'r'))
                config["RUN_TYPE"] = 1
                load_sphvars.write_sphvars(config,RUN_DIR+"/spinput.txt")

    print "SPAC run script ended"

if __name__ == "__main__":
    process_args()
    main()
