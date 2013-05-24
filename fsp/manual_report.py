#!/usr/bin/env python

""" 
Generates a latex report document designed to be added as an appendix to
the code manual. Designed to be run in the spac installation directory.

The full set of tests needs to have been run first (testall.py)

I gave up on letting Latex position everything and am using the textpos
package. For an A4 page is is simple just to make the grid 210x297
and position blocks to the nearest millimetre.

Positioning is {width}(x,y)?

Margins are 30mm each side and 50mm at the bottom.

There are three blocks at the moment:
Config text: {70}(30,30)

A raw string is defined as the document template, and we write the final
document by substituting for certain keywords in the template.

PIL is used to join the plots into a single image.

StringIO allows us to readlines() the template strings like a file, which is neat.

Andrew Charles
17 April 2009

"""

import sys
import os
from fsp_post import sph_results
from fsp_post import load_sphvars
import shutil
import glob
import netCDF4
import string
from PIL import Image
import StringIO

document_template = r""" 

\documentclass[10pt]{scrreprt}

\usepackage[absolute]{textpos}
\TPGrid[0mm,0mm]{210}{297}

\usepackage{graphicx}
\usepackage{wrapfig}
\usepackage[pdftex=true,colorlinks,bookmarks=true]{hyperref}

\begin{document}

<PAGE_ONE>

<PAGES>

\end{document}

"""


page_template = r"""

% see textpos manual - we need to \null to force tex to respect the newpage
\null\newpage

\begin{textblock}{80}(30,30)
<CONFIG>
\end{textblock}

\begin{textblock}{90}(110,30)
<SEQUENCE>
\end{textblock}

\begin{textblock}{90}(30,200)
<ENERGY>
\end{textblock}


"""


reportpath = 'test_cases/test_report'
ofname = 'test_report.tex'


def make_test_report(testpath='test_cases/test_output'):
    """ reportpath: where the report tex file and images live

    """

    if not os.path.isdir(testpath):
        print "Test directory does not exist, abort."
        sys.exit()

    if not os.path.isdir(reportpath):
        os.mkdir(reportpath)


    # We need a list of pages
    # Each page is a long string
    pages = []

    # We need a list of the cases to make a table of contents
    cases = []

    # Cycle through test directories
    # ------------------------------

    contents = os.listdir(testpath)
    for rundir in contents:
        if rundir[0] == '.':
            continue
        rundir = testpath + '/' + rundir
        if not os.path.isdir(rundir):
            continue
        runname = os.path.basename(rundir)
        print 'Generating report on ' + runname
        cases.append(runname)

        # Generate the config text
        # ------------------------
        nc_files = glob.glob(rundir + '/*.nc')
        print rundir
        print nc_files
        f = netCDF4.Dataset(nc_files[0],"r")

        config_list = ['Run Name:' + runname.replace('_','\_') + '\n\n'
                      ,r'\begin{description}' + '\n'
                      ,r'\setlength{\itemsep}{0pt}' + '\n'
                      ,r'\setlength{\parskip}{0pt}' + '\n'
                      ,r'\setlength{\parsep}{0pt}' + '\n' + '\n'
                      ]

        for att in f.ncattrs():
            config_list.append(r'\item [' + str(att).replace('_','\_') + '] ')
            config_list.append(str(getattr(f,str(att))).replace('_','\_') + '\n')

        config_list.append(r'\end{description}' + '\n')
        config_tex = string.join(config_list)

        # Generate the sequence text
        # --------------------------
        state_images = glob.glob(rundir + '/ra*.png')

        print 'Copying image files...'
        for img in state_images:
            dimgname = reportpath + '/' + runname + os.path.basename(img)
            shutil.copyfile(img,dimgname)

        # get the size
        im = Image.open(reportpath + '/' + runname + os.path.basename(state_images[0]))
        imx = im.size[0]
        imy = im.size[1]

        bimy = len(state_images) * im.size[1]
        imout = Image.new("RGB",(imx,bimy),(255,255,255))

        i = 0
        for img in state_images:
            dimgname = reportpath + '/' + runname + os.path.basename(img)
            im = Image.open(dimgname)
            frame = im.crop((0,0,im.size[0],im.size[1]))
            pco = ( (0,i*imy,imx,(i+1)*imy))
            imout.paste(frame,pco)
            i+=1

        imout.save(reportpath + '/sequence' + runname + '.png')

        sequence_list = [r'%\begin{wrapfigure}{r}{20mm}' + '\n'
                        ,r'\includegraphics[width=5cm]{sequence' + runname + '.png}\n'
                        ,r'%\end{wrapfigure}' + '\n']
        sequence_tex = string.join(sequence_list)

        # Generate the energy text
        # ------------------------

        # Energy conservation plot
        eimgname = runname + 'energy.png'
        shutil.copyfile(rundir + '/energy2.png',reportpath + '/' + eimgname)

        energy_list = [r'%\begin{wrapfigure}{r}{20mm}' + '\n'
                      ,r'\includegraphics[width=6cm]{' + eimgname + '}\n'
                      ,r'%\end{wrapfigure}' + '\n']

        energy_tex = string.join(energy_list)

        # Write out the page
        # ------------------

        page = []
        ptsio = StringIO.StringIO(page_template)


        for line in ptsio.readlines():
            line = line.replace("<CONFIG>",config_tex)
            line = line.replace("<SEQUENCE>",sequence_tex)
            line = line.replace("<ENERGY>",energy_tex)
            page.append(line)
            page_tex = string.join(page)

        # Add this page to the list of pages
        pages.append(page_tex)

    # Concatenate all the pages and substitute this into the document template
    pages_tex = string.join(pages)


    # Write out page 1
    # ----------------

    page_one = [r'\begin{itemize}'+'\n']
    for case in cases:
        page_one.append(r'\item ' + case.replace('_','\_')  + '\n')
    page_one.append(r'\end{itemize}'+'\n')
    page_one_tex = string.join(page_one)

    # Write the resulting document string to the output file
    # ------------------------------------------------------

    ofpath = reportpath + '/' + ofname
    ofile = open(ofpath,'w')

    for line in StringIO.StringIO(document_template).readlines():
        line = line.replace("<PAGE_ONE>",page_one_tex)
        line = line.replace("<PAGES>",pages_tex)
        ofile.write(line)
    ofile.close()


def make_pdf():
    cpath = os.path.abspath('.')
    os.chdir(reportpath)
    os.system("pdflatex " + ofname)
    os.chdir(cpath)


def main():
    """Main subroutine for smooth particle report writing"""
    make_test_report()
    make_pdf()

if __name__ == '__main__':
    main()
