def generate_html(data):
   #generates a simple html file
   ofile = open( 'output.html','w' )
   ofile.write( '<HTML><HEAD>')
   ofile.write( '<TITLE>sph variables</TITLE></HEAD>')
   ofile.write( '<BODY><H1>sph variables</H1>')
   ofile.write( '<p>')
   ofile.write( data )
   ofile.write( '</p>')
   ofile.write( '</BODY></HTML>')
