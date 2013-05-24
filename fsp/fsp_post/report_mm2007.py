#! /usr/bin/python
# This script moves through all subdirectories
# of the specified directory, and looks for
# sphvars.var files. It returns a report with the name
# contents.txt, containing selected information from
# sphvars.var
# Was written for the mm2007 conference, so it needs
# nice_plot images to be there
#            imgfname = 'nice_plot-sphstate.' + '%08d' %(largest) + '.png' 


import sys
import os
import sph_results
import load_sphvars
import shutil


#Check arguments
numargs = len(sys.argv)
print 'We have %d arguments' %numargs
print sys.argv[0]
if numargs>1:
   print sys.argv[1]
else:
   print "Expecting more arguments"
   sys.exit()

def main():
   """main subroutine contains high level flow"""
   print 'Looking for SP results in' + sys.argv[1]
   getContents(sys.argv[1])


def getContents(parent):
   argz = os.listdir(parent)
   ofile = open( 'contents.txt','w' )
   webfile = open( 'output.html','w' )
   webfile.write( '<HTML><HEAD>')
   webfile.write( '<TITLE>sph variables</TITLE></HEAD>')
   webfile.write( '<BODY><H1>sph variables</H1>')
   
   texfile = open( 'results.tex', 'w')
   
   #create a directory for the images
   if not(os.path.isdir('./images')):
      os.mkdir('./images')

   #for all folders and files in the specified path
   for a in argz:
      if os.path.isdir(parent + "/" + a):
         #for all directories in the path supplied...
         c=0
         if os.path.isfile(parent + "/" + a + "/sphvars.var"):
            #it's an sph results directory!
            
            spvf = open(parent + "/" + a + "/sphvars.var")
            params = load_sphvars.load_sphvars(spvf)
            tstep = load_sphvars.getTimestep(spvf)
            spvf.close()
      
            #get the id number of the run for multiple runs of the same temperature
            dirname = parent + "/" + a
            idnum = dirname[-1]

            #count the number of files in the directory
            dc = os.listdir(parent + "/" + a)
            for b in dc:
               c = c+1
               ostr = a + " %d %f" %(c,tstep)

            #print the output line to the terminal
            #and write it to file
            print ostr
            ostr = ostr+"\n"
            ofile.write( ostr )
            largest = sph_results.find_largest(parent +'/' + a)
            #find the maximum sphstate file
            print parent + "/" + a + ' %08d' %(largest)
               
            webfile.write( '<p>')      
            webfile.write(  'temp= %f' %(params['start_temp']) )
            webfile.write(  '<br/>run name= ' + a )
            webfile.write(  '<br/>steps=%d ' %(largest) )
            webfile.write(  '<br/>dt= %f' %(tstep) )
            webfile.write('</p>')
            
           #  nb this is start temp not thermostat temp. usually the same though
            
            #copy the final profile file to the images directory
            imgfname = 'nice_plot-sphstate.' + '%08d' %(largest) + '.png' 
            imgpath = parent +'/'+ a + '/' + imgfname
            print imgpath
            
            if os.path.isfile(imgpath):
               print 'is a file'
               imgarch = './images/' + '%d_' %(params['start_temp']*100) + idnum + '.png'
               shutil.copyfile(imgpath, imgarch)
               refpath = 'images/' + '%d_' %(params['start_temp']*100) + idnum + '.png'
               webfile.write('<img src="'+refpath+'" width="800" height="600"></img>')
               
              # texfile.write( '\\begin{figure}\n' )
             #  texfile.write( '\label{%d'%(params['start_temp']*100) + '}\n' )
               texfile.write( '\includegraphics[width=10cm]{'+refpath+'}\n' )
               texfile.write( '\\begin{verbatim}'+ refpath+ '\end{verbatim}\n' )
              # texfile.write( '\caption{%f' %(params['start_temp'])+ '}\n' )
             #  texfile.write( '\end{figure}\n\n' )

            
         else:
            # it's not an sph run directory
            print "non sph run object"

   ofile.close()
   webfile.write( '</BODY></HTML>')
   webfile.close()
   texfile.close()

#run main
main()
