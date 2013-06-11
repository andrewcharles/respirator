#! /usr/bin/python
"""  mpegify.py 
 This takes directories with insane amounts of pngs in them
 and makes movies, at more or less a set frame rate, of images with the
 same starting string

  TODO: pass in the glob string


"""
import glob
import os

MENCODER=0
FF=1
 
def animate(globstring,fps):
    """makes an mpeg"""
    print 
    if MENCODER:
        cmdstring = "mencoder 'mf://" + globstring + "*.png' " \
                    + "-nosound " \
                    + "-mf type=png:fps=%d " %fps  \
                    + "-ovc lavc " \
                    + "-lavcopts vcodec=msmpeg4v2:vbitrate=1200:cmp=2:mbd=2 "\
                    + "-oac copy -o " \
                    + globstring + ".avi"

    elif FF:
        cmdstring = "ffmpeg -y -qscale 3 -i " + globstring + "%05dI.png " + \
                    globstring + ".mpg"

    print cmdstring
    os.system(cmdstring)

def process(forms):
    for f in forms:
        length = len(glob.glob(f+"*.png"))
        if length < 5:
            fps =5 
        elif length < 20:
            fps = 5
        else:
            fps = 10
        #make the movies
        animate(f,fps)

# Main routine here!
#-------------------
def old_main():
    # 1. Take stock of what pngs we have here
    all_pngs = glob.glob('*.png')
    print "there are " + str(len(all_pngs)) + " pngs in here"
    # use a python set to get only unique first five letters?
    forms = ["a","d","t","p","r"] 
        #os.remove(all_pngs)
    process(forms)

if __name__ == "__main__":
    #main()
    animate('f',12)
