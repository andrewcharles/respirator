"""
    Viewer module for particle system
    Copyright Andrew Charles 2008
    All rights reserved.

    Requires Pyglet.

    To mess around with the various visualisation
    options edit the function redraw()

    The trackball camera is actually a bit non-intuitive for me.
    I think I would prefer more fpsy controls:

    W,S - zoom in,out
    A,D - pan left,right
    Q,E - rotate left,right


"""

import sys
import particles
import pyglet
pyglet.options['debug_gl'] = True
from pyglet.gl import *
import numpy
from pygarrayimage.arrayimage import ArrayInterfaceImage
from pylab import *
import pdb
import spkernel
import interpolate
import renderspam
import scipy
import scipy.ndimage
import pylab
from time import time
#from PIL import Image as pilmage
#import scipy.misc.pilutil
#from trackball_camera import TrackballCamera


class ParticleInfo:
    """ Display particle system information."""
    """ The motivation for this class is to make information like
        the time to execute certain subroutines available to a
        range of objects, and keep the pyglet window for the
        rendered view.
    """

    def __init__(self,p):
        self.win = pyglet.window.Window(200,300
                 ,visible=True
                 ,resizable=True
                 ,caption='Pyticles info')
        self.labels = []
        self.npartlab = pyglet.text.Label("np label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=10,y=440 )
        self.labels.append(self.npartlab)

#class ViewController:
    #""" Pyglet key mappings for controlling the camera and other
        #elements of the view.
    #"""

class plabel(pyglet.text.Label):
    """ Subclass of pyglet label with some defaults selected."""
    def __init__(self,name,loc=(0,0),col=(255,0,0,255)):
        pyglet.text.Label.__init__(self,
            name,
            font_name="Arial",
            font_size=12,
            color=col,x=loc[0],y=loc[1] )


WINDOW_WIDTH = 640
WINDOW_HEIGHT = 480
RES = 1.0

class ParticleView2d:
    " A particle viewer"
    def __init__(self,particles):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT
                 ,visible=False,caption='Pyticles')
        self.fps = 0
        self.fpslab = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=4 )
        self.npart = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=4 )
        self.nebs = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=18 )
        self.maxvol = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=18 )
        
        self.win.set_visible()
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / float(particles.box.xmax)
        self.ymap = WINDOW_HEIGHT / float(particles.box.ymax)
        self.fps = 0

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
        radius=5
        for i in range(p.n):
            r = p.r[i,0] * self.xmap, p.r[i,1] * self.ymap
            glBegin(GL_POLYGON)
            glColor3f(p.colour[0],p.colour[1],p.colour[2])
            for angle in range(6):
                a = radians(angle*60)
                glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            glEnd()

    def draw_neighbours(self,p):
        """ issues the opengl commands to draw 
            lines connecting neighbours """
        for i in range (p.nl_default.nip):
            glBegin(GL_LINES)
            glColor3f(0.3,0.3,0.0)
            a = p.nl_default.iap[i,0]
            b = p.nl_default.iap[i,1]
            glVertex2f(p.r[a,0]*self.xmap,p.r[a,1]*self.ymap)
            glVertex2f(p.r[b,0]*self.xmap,p.r[b,1]*self.ymap)
            glEnd()
  
    def hud(self,p):
        """ Draws the heads up display """
        k = 20
        xi = 10
        yi = 450
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.fpslab.draw()
        self.nebs.text = "n: "+str(p.n)
        self.nebs.draw()
        self.npart.text = "pairs: "+str(p.nl_default.nip)
        self.npart.draw()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        self.hud(p)
        self.draw_neighbours(p)
        self.draw_particles(p)


class SmoothParticleView2d(ParticleView2d):
    " A particle viewer"
    def __init__(self):
        #config = Config(alpha_size=8)
        self.win = pyglet.window.Window(WINDOW_WIDTH,WINDOW_HEIGHT,visible=False,caption='Pyticles')
        self.fps = 0
        self.fpslab = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=4 )
        self.npart = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=4 )
        self.nebs = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=10,y=18 )
        self.maxvol = pyglet.text.Label("Hello pyglet",font_name="Arial", \
        font_size=12,color =(255,0,0,255),x=100,y=18 )
        
        self.win.set_visible()
        # multipliers to map system to screen coordinates
        self.xmap = WINDOW_WIDTH / float(particles.XMAX)
        self.ymap = WINDOW_HEIGHT / float(particles.YMAX)
        self.img = pyglet.image.load('p.png')
        self.fps = 0
        # gridx,gridy is a grid at the rendering resolution
        self.gridx,self.gridy = renderspam.get_grid_map(0.0,particles.XMAX, \
        0.0,particles.YMAX,RES)
        # rx,ry is a pixel grid
        rx,ry = scipy.mgrid[0:WINDOW_WIDTH:1,0:WINDOW_HEIGHT:1]
        dx = RES*WINDOW_WIDTH/float(particles.XMAX)
        dy = RES*WINDOW_HEIGHT/float(particles.YMAX)
        #xvals = (rx-dx/RES)/dx #(rx-RES?)
        xvals = (rx-dx/2)/dx #(rx-RES?)
        yvals = (ry-dy/2)/dy
        self.G = numpy.array([xvals,yvals])
          

    def render_density(self,p):
        """ Generates a splashmap, then creates an image the
            size of the viewing window and interpolates the splashmap
            onto it.
            Currently dependent on vasp modules. 
        """ 
        glColor3f(1.0,1.0,1.0)
        #cutoff=p.nl_default.cutoff_radius
        cutoff=5
        bounds = 0,particles.XMAX,0,particles.YMAX
        Z = interpolate.splash_map(self.gridx,self.gridy,p.r[0:p.n,:], \
            p.m[0:p.n],p.rho[0:p.n],p.rho[0:p.n],p.h[0:p.n],bounds,cutoff)
        D = scipy.ndimage.map_coordinates(Z,self.G,order=1,mode='nearest', \
            prefilter=False)
        #D = numpy.random.random([640,480])
        # need to create an image the size of the screen
        # and then interpolate based on our splashmap
        D = numpy.rot90(D)
        #print amax(D/amax(D))
        D = 255.*(D/amax(D))
        D = numpy.flipud(D)
        D = numpy.cast['uint8'](D)
        A = ArrayInterfaceImage(D)
        img = A.image_data
        #A = pilmage.fromarray(D,'L')
        #img = pyglet.image.load('wtf.png',A)
        img.blit(0,0)

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
        radius=5
        for i in range(p.n):
        #    self.img.blit(p.r[i,0],p.r[i,1])
            r = p.r[i,0] * self.xmap, p.r[i,1] * self.ymap
            glBegin(GL_POLYGON)
            glColor3f(p.colour[0],p.colour[1],p.colour[2])
            for angle in range(6):
                a = radians(angle*60)
                glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            glEnd()

    def draw_particle_sprites(self,p):
        """ Draws particles as sprites.
        """
        for i in range(p.n):
            self.img.blit(self.xmap*p.r[i,0],self.ymap*p.r[i,1])

        
    def draw_neighbours(self,p):
        """ issues the opengl commands to draw 
            lines connecting neighbours """
        for i in range (p.nl_default.nip):
            glBegin(GL_LINES)
            glColor3f(0.3,0.3,0.0)
            a = p.nl_default.iap[i,0]
            b = p.nl_default.iap[i,1]
            glVertex2f(p.r[a,0]*self.xmap,p.r[a,1]*self.ymap)
            glVertex2f(p.r[b,0]*self.xmap,p.r[b,1]*self.ymap)
            glEnd()
  
    def hud(self,p):
        """ Draws the heads up display """
        k = 20
        xi = 10
        yi = 450
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.fpslab.draw()
        self.nebs.text = "n: "+str(p.n)
        self.nebs.draw()
        self.npart.text = "pairs: "+str(p.nl_default.nip)
        self.npart.draw()
        self.maxvol.text = "max vol: %4.2f" %(max(p.m[0:p.n]/p.rho[0:p.n]))
        self.maxvol.draw()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        #self.win.dispatch_events()
        #glClear(GL_COLOR_BUFFER_BIT)
        #glLoadIdentity()
        #self.render_density(p)
        self.hud(p)
        #self.draw_neighbours(p)
        self.draw_particles(p)
        #self.draw_particle_sprites(p)



class ParticleView:
    " A particle viewer"
    def __init__(self,p,wsize=(640,480,480),box=(250,250,250)):
        #config = Config(alpha_size=8)

        self.WINDOW_WIDTH = wsize[0]
        self.WINDOW_HEIGHT = wsize[1]
        self.WINDOW_DEPTH = wsize[2]

        # These are the opengl box dimensions
        # We need a scaling factor to scale the simulation
        # box dimension to opengl coordinates.
        # self.xmap =  box_width / particles.xmap
        self.BOX_WIDTH = box[0]
        self.BOX_HEIGHT = box[1]
        self.BOX_DEPTH = box[2]
        self.PSIZE = 10.0
        self.RES = 1.1

        self.win = pyglet.window.Window(self.WINDOW_WIDTH,self.WINDOW_HEIGHT
                 ,visible=False,caption='Pyticles')

        self.quadric = gluNewQuadric()
        gluQuadricDrawStyle(self.quadric,GLU_FILL)
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_CULL_FACE)
        
        self.bg_color = (0.8,0.8,0.8,1.0)

        # Set up labels. These labels are pretty tightly coupled to the
        # gui buttons
        self.fps = 0

        self.labels=[]

        self.fpslab = plabel("fps label",loc=(10,460))
        self.npartlab = plabel("np label",loc=(10,440))
        self.nebslab = plabel("ni label",loc=(10,420))
        self.dtlab = plabel("dt label",loc=(10,400))
        self.eyelab = plabel("eye label",loc=(10,380) )
        self.drawtimelab = plabel("drawtime label", loc=(10,360) )
        self.derivtimelab = plabel("derivtime label",loc=(10,340) )
        self.integtimelab = plabel("integtime label",loc=(10,320) )
        self.stepslab = plabel("steps label",loc=(10,300) )
        self.updatetimelab = plabel("updatetime label",loc=(10,280) )

        self.pairseptimelab = plabel("pairseptime label",loc=(10,260)
            ,col=(255,0,255,255) )
        self.spamtimelab = plabel("spamtime label",loc=(10,240)
            ,col=(255,0,255,255) )
        self.forcetimelab = plabel("forcetime label",loc=(10,220)
            ,col=(255,0,255,255) )
        self.energylab = plabel("kinetic energy",loc=(10,200)
            ,col=(255,0,255,255) )
        
        self.labels.append(self.fpslab)
        self.labels.append(self.npartlab)
        self.labels.append(self.nebslab)
        self.labels.append(self.dtlab)
        self.labels.append(self.eyelab)
        self.labels.append(self.drawtimelab)
        self.labels.append(self.forcetimelab)
        self.labels.append(self.derivtimelab)
        self.labels.append(self.pairseptimelab)
        self.labels.append(self.integtimelab)
        self.labels.append(self.stepslab)
        self.labels.append(self.updatetimelab)
        self.labels.append(self.spamtimelab)
        self.labels.append(self.energylab)

        self.tbcam  = None #TrackballCamera(200.0) 
        self.win.on_resize = self.resize
        self.win.on_mouse_press = self.on_mouse_press
        self.win.on_mouse_drag = self.on_mouse_drag
        self.win.on_key_press = self.wsad_press
        
        self.win.set_visible()
        
        # multipliers to map system to opengl box coordinates
        self.xmap = self.BOX_WIDTH / float(p.box.xmax)
        self.ymap = self.BOX_HEIGHT / float(p.box.ymax)
        self.zmap = self.BOX_DEPTH / float(p.box.zmax)

        # Where is the center of the rectangle?
        cx = float(particles.XMAX/2.0) * self.xmap
        cy = float(particles.YMAX/2.0) * self.ymap
        cz = float(particles.ZMAX/2.0) * self.zmap
        self.tbcam.cam_focus = (cx,cy,cz)

        # zoom out
        self.tbcam.cam_eye[0] = 490
        self.tbcam.cam_eye[1] = 430
        self.tbcam.cam_eye[2] = 400
        self.tbcam.update_modelview()

        self.zoom = 0
        self.rotation = 0

        # hit the movement buttons to start rotations etc
        self.eyespeed = [0.0,0.0,0.0]

        """ Variables for measuring performance. """
        self.timing = {}
        self.timing['Draw time'] = -1

    def center_eye(self):
        self.tbcam.cam_eye[2] = -600
        self.tbcam.cam_eye[1] = -0
        self.tbcam.cam_eye[0] = -0
        self.tbcam.update_modelview()
        self.eyespeed = [0.0,0.0,0.0]

    def draw_box(self):
        x = self.BOX_WIDTH
        y = self.BOX_HEIGHT
        z = self.BOX_DEPTH

        """
    
    (0,y,z) --------- (x,y,z)
       |\                | \
       | \          3    |  \
       |  \              |   \
       |   \             |    \
       |   (0,y,0) -------- (x,y,0)
       | 1  |            |      |
    (0,0,z) | ------- (x,0,z)   |
        \   |       2       \   |
         \  |                \  |
          \ |                 \ |
           \|                  \|
        (0,0,0) ----------- (x,0,0)

        """
        
        glBegin(GL_LINE_STRIP)
        glColor3f(0.0,0.0,1.0)
       
        # Face 1 
        glVertex3f(0,0,0)
        glVertex3f(0,0,z)
        glVertex3f(0,y,z)
        glVertex3f(0,y,0)
        glVertex3f(0,0,0)

        # Face 2
        glVertex3f(x,0,0)
        glVertex3f(x,y,0)
        glVertex3f(0,y,0)

        # Face 3
        glVertex3f(0,y,z)
        glVertex3f(x,y,z)
        glVertex3f(x,y,0)
        glVertex3f(0,y,0)

        # Face 5
        glVertex3f(0,y,z)
        glVertex3f(0,0,z)
        glVertex3f(x,0,z)
        glVertex3f(x,y,z)
        glVertex3f(0,y,z)

        # Face 6
        glVertex3f(x,y,z)
        glVertex3f(x,y,0)
        glVertex3f(x,0,0)
        glVertex3f(x,0,z)

        glEnd()
        
        
    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        radius=5
        for i in range(p.n):
            r = p.r[i,0] * self.xmap \
              , p.r[i,1] * self.ymap \
              , p.r[i,2] * self.zmap
            glColor3f(1.0, 1.0/(i+1), 0)
            glPushMatrix()
            glTranslatef(r[0],r[1],r[2])
            gluSphere(self.quadric,self.PSIZE,4,4)
            
            #glBegin(GL_POLYGON)
            #glColor3f(p.colour[0],p.colour[1],p.colour[2])
            #for angle in range(4):
            #    a = radians(angle*60)
            #    glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            #glEnd()
            
            glPopMatrix()


    def draw_neighbours(self,p):
        """ Issues the opengl commands to draw 
            lines connecting neighbours. """
        for i in range (p.nl_default.nip):
            glBegin(GL_LINES)
            glColor3f(0.3,0.3,0.3)
            a = p.nl_default.iap[i,0]
            b = p.nl_default.iap[i,1]
            glVertex3f(
                p.r[a,0]*self.xmap,
                p.r[a,1]*self.ymap,
                p.r[a,2]*self.zmap)
            glVertex3f(
                p.r[b,0]*self.xmap,
                p.r[b,1]*self.ymap,
                p.r[b,2]*self.zmap)
            glEnd()
  
    def hud(self,p):
        """ Draws the heads up display """
        k = 20
        xi = 10
        yi = 450

        self.set_ortho()
        glPushMatrix()
        glLoadIdentity()
       
        eye = self.tbcam.cam_eye
        self.eyelab.text = 'eye: %d, %d, %d' %(eye[0],eye[1],eye[2])
        self.dtlab.text = "dt: %4.2f" %(p.dt)
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.nebslab.text = "n: %3d" %(p.n)
        nebs = 0
        for nl in p.nlists:
            nebs += nl.nip
        self.npartlab.text = "pairs: %3d" %(nebs)
        self.drawtimelab.text = "drawtime:  %5.3f" %(self.timing['Draw time'])
        self.forcetimelab.text = \
            "forcetime:   %7.5f" %(p.timing['force time'])
        self.derivtimelab.text = "derivtime: %5.3f" %(p.timing['deriv time'])
        self.pairseptimelab.text = \
            "pairtime:     %7.5f" %(p.timing['pairsep time'])
        self.integtimelab.text = \
            "integtime: %5.3f" %(p.timing['integrate time'])
        self.spamtimelab.text = "spamtime:  %7.5f" %(p.timing['SPAM time'])
        self.stepslab.text = "steps:     %5d" %(p.steps)
        self.updatetimelab.text = \
            "updatetime: %5.3f" %(p.timing['update time'])
        self.energylab.text = \
            "EK: %5.3f" %( np.mean(0.5*(np.linalg.norm(p.v)**2)*p.m) )

        for lab in self.labels:
            lab.draw()
        glFlush()
        glPopMatrix()

        self.reset_perspective()

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glClearColor(*self.bg_color)
        #glLoadIdentity()
        self.win.clear()
        
    def redraw(self,p):
        t = time()
        self.win.dispatch_events()
        glClearColor(*self.bg_color)
        #self.win.clear()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #self.draw_neighbours(p)
        self.draw_particles(p)
        self.draw_box()
        self.hud(p)
        self.timing['Draw time'] = time() - t

    def set_ortho(self):
        """ Sets up an orthographic projection for drawing text and the like.
        """
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        #gluOrtho2D(0,width,0,height)
        glOrtho(0,self.WINDOW_WIDTH,0,self.WINDOW_HEIGHT,-1,1)
        glMatrixMode(GL_MODELVIEW)

    def reset_perspective(self):
        """ Pops the projection matrix. Used for getting the
            3D perspective view back.
        """
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        self.tbcam.update_modelview()

    def resize(self,width, height):
        """Setup 3D projection for window"""
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, 1.0 * width/height, 0.001, 10000)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.tbcam.update_modelview()

    def norm1(self,x,maxx):
        """given x within [0,maxx], scale to a range [-1,1]."""
        return (2.0 * x - float(maxx)) / float(maxx)

    def on_mouse_press(self, x, y, button, modifiers):
        if button == pyglet.window.mouse.LEFT:
            self.tbcam.mouse_roll(
                self.norm1(x,self.WINDOW_WIDTH),
                self.norm1(y,self.WINDOW_HEIGHT),
                False)
        elif button == pyglet.window.mouse.RIGHT:
            self.tbcam.mouse_zoom(
                self.norm1(x,self.WINDOW_WIDTH),
                self.norm1(y,self.WINDOW_HEIGHT),
                False)

    def on_mouse_drag(self,x, y, dx, dy, buttons, modifiers):
        if buttons & pyglet.window.mouse.LEFT:
            if modifiers & pyglet.window.key.MOD_SHIFT:
                self.tbcam.mouse_roll(
                    self.norm1(x,self.WINDOW_WIDTH),
                    self.norm1(y,self.WINDOW_HEIGHT))
            else:
                self.tbcam.mouse_zoom(
                    self.norm1(x,self.WINDOW_WIDTH),
                    self.norm1(y,self.WINDOW_HEIGHT))

    def update_eye(self,t):
            self.tbcam.cam_eye[0] += self.eyespeed[0]
            self.tbcam.cam_eye[1] += self.eyespeed[1]
            self.tbcam.cam_eye[2] += self.eyespeed[2]
            self.tbcam.update_modelview()

    def wsad_press(self,symbol,modifiers):
        """ zoom and rotate directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
            self.tbcam.mouse_zoom_diff(0.2,0.2)
            # zoom in
            #self.eyespeed[2] += 100
        if symbol == pyglet.window.key.S:
            self.tbcam.mouse_zoom_diff(-0.2,-0.2)
            #self.eyespeed[2] -= 100
        if symbol == pyglet.window.key.A:
            self.eyespeed[1] -= 30
        if symbol == pyglet.window.key.D:
            self.eyespeed[1] += 30
        if symbol == pyglet.window.key.Q:
            self.eyespeed[0] -= 30
        if symbol == pyglet.window.key.E:
            self.eyespeed[0] += 30
        if symbol == pyglet.window.key.X:
            self.center_eye()
        if symbol == pyglet.window.key.C:
            self.eyespeed[0] = 0
            self.eyespeed[1] = 0
            self.eyespeed[2] = 0

    def wsad_release(self,symbol,modifiers):
        """ zoom, pan and rotate directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
            pass
        if symbol == pyglet.window.key.S:
            pass
        if symbol == pyglet.window.key.A:
            pass
        if symbol == pyglet.window.key.D:
            pass
            

class SmoothParticleView(ParticleView):
    """ Extends ParticleView to provide specialised visualisation of
        smooth particle systems. """
    def __init__(self,p,wsize=(640,480,480),box=(250,250,250)):
        #ParticleView.__init__(self,p,wsize=(640,480,480),box=(250,250,250))
        ParticleView.__init__(self,p)
        self.WINDOW_WIDTH = wsize[0]
        self.WINDOW_HEIGHT = wsize[1]
        self.WINDOW_DEPTH = wsize[2]

        self.BOX_WIDTH = box[0]
        self.BOX_HEIGHT = box[1]
        self.BOX_DEPTH = box[2]
        self.PSIZE = 10.0
        self.RES = 1.1

        self.maxvol = pyglet.text.Label("maxvol label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=10,y=100 )

        # Set up for two dimensional density rendering
        # gridx,gridy is a grid at the rendering resolution
        # rx,ry is a pixel grid
        self.gridx,self.gridy = renderspam.get_grid_map(0.0,p.box.xmax, \
        0.0,p.box.ymax,self.RES)
        rx,ry = scipy.mgrid[0:self.BOX_WIDTH:1,0:self.BOX_HEIGHT:1]
        dx = self.RES*self.BOX_WIDTH/float(p.box.xmax)
        dy = self.RES*self.BOX_HEIGHT/float(p.box.ymax)
        #xvals = (rx-dx/RES)/dx #(rx-RES?)
        xvals = (rx-dx/2)/dx #(rx-RES?)
        yvals = (ry-dy/2)/dy
        self.G = numpy.array([xvals,yvals])
         
    def render_density(self,p):
        """ Generates a splashmap, then creates an image the
            size of the viewing window and interpolates the splashmap
            onto it.
            Currently dependent on vasp modules.
            Written as it is for two dimensions.
        """
        glColor3f(1.0,1.0,1.0)
        #cutoff=p.nl_default.cutoff_radius
        cutoff=5
        bounds = 0,p.box.xmax,0,p.box.ymax
        Z = interpolate.splash_map_3d(
            self.gridx,
            self.gridy,
            0.0,
            p.r[0:p.n,:], 
            p.m[0:p.n],
            p.rho[0:p.n],
            p.rho[0:p.n],
            p.h[0:p.n],
            bounds,cutoff)

        D = scipy.ndimage.map_coordinates(Z,self.G,order=1,mode='nearest', \
            prefilter=False)
        #D = numpy.random.random([640,480])
        # need to create an image the size of the screen
        # and then interpolate based on our splashmap
        D = numpy.rot90(D)
        #print amax(D/amax(D))
        if amax(D) == 0.0:
            D = D
        else:
            D = 255.*(D/amax(D))
        D = numpy.flipud(D)
        D = numpy.cast['uint8'](D)
        A = ArrayInterfaceImage(D)
        img = A.image_data
        #A = pilmage.fromarray(D,'L')
        #img = pyglet.image.load('wtf.png',A)
        img.blit(0,0)

    def draw_particle_sprites(self,p):
        """ Draws particles as sprites.
        """
        for i in range(p.n):
            self.img.blit(self.xmap*p.r[i,0],self.ymap*p.r[i,1])

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        radius=5
        for i in range(p.n):
            r = p.r[i,0] * self.xmap \
              , p.r[i,1] * self.ymap \
              , p.r[i,2] * self.zmap
            a = p.t[i]/2.0
            glColor3f(0, 0, a)
            glPushMatrix()
            glTranslatef(r[0],r[1],r[2])
            #gluSphere(self.quadric,self.PSIZE,10,4)
            
            glBegin(GL_POLYGON)
            glColor3f(p.colour[0],p.colour[1],p.colour[2])
            for angle in range(4):
                a = radians(angle*60)
                glVertex2f(r[0]+sin(a)*radius,r[1]+cos(a)*radius)
            glEnd()
            
            glPopMatrix()

    def redraw(self,p):
        t = time()
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #self.render_density(p)
        #self.hud(p)
        #self.draw_neighbours(p)
        #self.draw_particles(p)
        self.draw_box()
        self.timing['Draw time'] = time() - t


class ZPRView(SmoothParticleView):
    """ Development class for ortho projection and zpr instead of
        trackball.
    """
    def __init__(self,p,wsize=(640,480,480),box=(250,250,250),lighting=True):
        #ParticleView.__init__(self,p)
        self.WINDOW_WIDTH = wsize[0]
        self.WINDOW_HEIGHT = wsize[1]
        self.WINDOW_DEPTH = wsize[2]

        self.BOX_WIDTH = box[0]
        self.BOX_HEIGHT = box[1]
        self.BOX_DEPTH = box[2]
        self.PSIZE = 10.0
        self.RES = 1.1

        self.win = pyglet.window.Window(self.WINDOW_WIDTH,self.WINDOW_HEIGHT
                 ,visible=False,caption='Pyticles')

        self.sphere = gluNewQuadric()
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_CULL_FACE)
        col = [1.0,1.0,1.0]
        color = (GLfloat * 3)(*col)

        if lighting:
            glEnable(GL_LIGHTING)
            glEnable(GL_LIGHT0)
            glLightfv(GL_LIGHT0, GL_AMBIENT,color)
            gluQuadricDrawStyle(self.sphere,GLU_FILL)
        else:
            gluQuadricDrawStyle(self.sphere,GLU_LINE)
        
        self.bg_color = (0.8,0.8,0.8,1.0)

        self.maxvol = pyglet.text.Label("maxvol label",font_name="Arial", \
            font_size=12,color =(255,0,0,255),x=10,y=100 )

        # Set up for two dimensional density rendering
        # gridx,gridy is a grid at the rendering resolution
        # rx,ry is a pixel grid
        self.gridx,self.gridy = renderspam.get_grid_map(0.0,p.box.xmax, \
        0.0,p.box.ymax,self.RES)
        rx,ry = scipy.mgrid[0:self.BOX_WIDTH:1,0:self.BOX_HEIGHT:1]
        dx = self.RES*self.BOX_WIDTH/float(p.box.xmax)
        dy = self.RES*self.BOX_HEIGHT/float(p.box.ymax)
        #xvals = (rx-dx/RES)/dx #(rx-RES?)
        xvals = (rx-dx/2)/dx #(rx-RES?)
        yvals = (ry-dy/2)/dy
        self.G = numpy.array([xvals,yvals])

        # Set up labels. These labels are pretty tightly coupled to the
        # gui buttons
        self.fps = 0

        self.labels=[]

        #for lname in ['fps label','np label','ni label','dt label','eye label',
        #    'drawtime label','forcetime label','derivtime label',
        #    'pairseptime label']

        self.fpslab = plabel("fps label",loc=(10,460))
        self.npartlab = plabel("np label",loc=(10,440))
        self.nebslab = plabel("ni label",loc=(10,420))
        self.dtlab = plabel("dt label",loc=(10,400))
        self.eyelab = plabel("eye label",loc=(10,380) )
        self.drawtimelab = plabel("drawtime label", loc=(10,360) )
        self.derivtimelab = plabel("derivtime label",loc=(10,340) )
        self.integtimelab = plabel("integtime label",loc=(10,320) )
        self.stepslab = plabel("steps label",loc=(10,300) )
        self.updatetimelab = plabel("updatetime label",loc=(10,280) )

        self.pairseptimelab = plabel("pairseptime label",loc=(10,260)
            ,col=(255,0,255,255) )
        self.spamtimelab = plabel("spamtime label",loc=(10,240)
            ,col=(255,0,255,255) )
        self.forcetimelab = plabel("forcetime label",loc=(10,220)
            ,col=(255,0,255,255) )
        self.templab = plabel("energy label",loc=(10,200)
            ,col=(255,0,255,255) )
        self.energylab = plabel("energy label",loc=(10,180)
            ,col=(255,0,255,255) )
        self.densitylab = plabel("density label",loc=(10,160)
            ,col=(255,0,255,255) )

        self.labels.append(self.fpslab)
        self.labels.append(self.npartlab)
        self.labels.append(self.nebslab)
        self.labels.append(self.dtlab)
        self.labels.append(self.eyelab)
        self.labels.append(self.drawtimelab)
        self.labels.append(self.forcetimelab)
        self.labels.append(self.derivtimelab)
        self.labels.append(self.pairseptimelab)
        self.labels.append(self.integtimelab)
        self.labels.append(self.stepslab)
        self.labels.append(self.updatetimelab)
        self.labels.append(self.spamtimelab)
        self.labels.append(self.templab)
        self.labels.append(self.energylab)
        self.labels.append(self.densitylab)

        self.win.on_resize = self.resize
        self.win.on_mouse_press = self.on_mouse_press
        self.win.on_mouse_drag = self.on_mouse_drag
        self.win.on_key_press = self.wsad_press
        
        self.win.set_visible()
        
        # multipliers to map system to opengl box coordinates
        self.xmap = self.BOX_WIDTH / float(p.box.xmax)
        self.ymap = self.BOX_HEIGHT / float(p.box.ymax)
        self.zmap = self.BOX_DEPTH / float(p.box.zmax)

        # Where is the center of the rectangle?
        self.cx = float(particles.XMAX/2.0) * self.xmap
        self.cy = float(particles.YMAX/2.0) * self.ymap
        self.cz = float(particles.ZMAX/2.0) * self.zmap

        self.zoom = 0
        self.rotation = 0

        """ Variables for measuring performance. """
        self.timing = {}
        self.timing['Draw time'] = -1

    def center_eye(self):
        pass

    def hud(self,p):
        """ Draws the heads up display """
        k = 20
        xi = 10
        yi = 450

        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0,self.WINDOW_WIDTH,0,self.WINDOW_HEIGHT,0,10)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
       
        self.eyelab.text = 'Eye coordinates.'
        self.dtlab.text = "dt: %4.2f" %(p.dt)
        self.fpslab.text = "fps: %4.2f" %(self.fps)
        self.nebslab.text = "n: %3d" %(p.n)
        nebs = 0
        for nl in p.nlists:
            nebs += nl.nip
        self.npartlab.text = "pairs: %3d" %(nebs)
        self.drawtimelab.text = \
            "drawtime:  %5.3f" %(self.timing['Draw time'])
        self.forcetimelab.text = \
            "forcetime:   %7.5f" %(p.timing['force time'])
        self.derivtimelab.text = "derivtime: %5.3f" %(p.timing['deriv time'])
        self.pairseptimelab.text = \
            "pairtime:     %7.5f" %(p.timing['pairsep time'])
        self.integtimelab.text =   \
            "integtime: %5.3f" %(p.timing['integrate time'])
        self.spamtimelab.text = "spamtime:  %7.5f" %(p.timing['SPAM time'])
        self.stepslab.text = "steps:     %5d" %(p.steps)
        self.templab.text = "mean temp:     %5.3f" %(np.mean(p.t))
        self.updatetimelab.text = \
            "updatetime: %5.3f" %(p.timing['update time'])
        self.energylab.text = \
            "EK: %5.3f" %( np.mean(0.5*(np.linalg.norm(p.v)**2)*p.m) )
        self.densitylab.text = "rhomax: %5.3f" %(p.rho.max())
        for lab in self.labels:
            lab.draw()
        glFlush()
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

    def clear(self): 
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT)
        glClearColor(*self.bg_color)
        self.win.clear()
        
    def redraw(self,p):
        t = time()
        self.win.dispatch_events()
        glClearColor(*self.bg_color)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.draw_particles(p)
        self.draw_box()
        self.hud(p)
        self.timing['Draw time'] = time() - t

    def resize(self,width, height):
        """Setup 3D projection for window"""
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        left = -10
        right = self.BOX_WIDTH +10
        top = self.BOX_HEIGHT + 10
        bottom = -10
        znear = -1000
        zfar = 10000
        glOrtho(left,right,bottom,top,znear,zfar)
        glMatrixMode(GL_MODELVIEW)
        #glLoadIdentity()

    def norm1(self,x,maxx):
        """given x within [0,maxx], scale to a range [-1,1]."""
        return (2.0 * x - float(maxx)) / float(maxx)

    def on_mouse_press(self, x, y, button, modifiers):
        if button == pyglet.window.mouse.LEFT:
            print 'rotate'
        elif button == pyglet.window.mouse.RIGHT:
            print 'zoom'

    def on_mouse_drag(self,x, y, dx, dy, buttons, modifiers):
        if buttons & pyglet.window.mouse.LEFT:
            if modifiers & pyglet.window.key.MOD_SHIFT:
                print 'shift drag'
            else:
                print 'zoom'

    def wsad_press(self,symbol,modifiers):
        """ zoom and rotate directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
           glTranslatef(self.cx,self.cy,self.cz)
           glScalef(1.1,1.1,1.1)
           glTranslatef(-self.cx,-self.cy,-self.cz)
        if symbol == pyglet.window.key.S:
           glTranslatef(self.cx,self.cy,self.cz)
           glScalef(0.9,0.9,0.9)
           glTranslatef(-self.cx,-self.cy,-self.cz)
        if symbol == pyglet.window.key.A:
            glRotatef(10,0,1.0,0) 
        if symbol == pyglet.window.key.D:
            glRotatef(-10,0.0,1.0,0) 
        if symbol == pyglet.window.key.Q:
            glRotatef(5,1.0,0.0,0) 
        if symbol == pyglet.window.key.E:
            glRotatef(-5,1.0,0.0,0) 
        if symbol == pyglet.window.key.X:
            print 'zoom out'
        if symbol == pyglet.window.key.C:
            print 'zoom out'

    def wsad_release(self,symbol,modifiers):
        """ zoom, pan and rotate directly with the fps keys.
        """
        if symbol == pyglet.window.key.W:
            pass
        if symbol == pyglet.window.key.S:
            pass
        if symbol == pyglet.window.key.A:
            pass
        if symbol == pyglet.window.key.D:
            pass

    def draw_particles(self,p):
        """ issues the opengl commands to draw the 
            particle system """
        radius=5
        scale = 1.0
        for i in range(p.n):
            r = p.r[i,0] * self.xmap \
              , p.r[i,1] * self.ymap \
              , p.r[i,2] * self.zmap
            a = p.rho[i]/2.0
            glColor3f(0, 0, a)
            glPushMatrix()
            glTranslatef(r[0],r[1],r[2])
            gluSphere(self.sphere,self.PSIZE,10,4)
            glBegin(GL_LINES);
            glColor3f(1.0,0.0,1.0);
            #glVertex3f(r[0],r[1],r[2]);
            glVertex3f(0,0,0);
            a = p.vdot[i,0] * self.xmap * scale\
              , p.vdot[i,1] * self.ymap * scale\
              , p.vdot[i,2] * self.zmap * scale
            glVertex3f(a[0],a[1],a[2])
            glEnd()
            glPopMatrix()

    def density_histogram(self,p):
        """ Draws the heads up display """
        hx = (0,200)
        hy = (0,100)
        hdx = 10.
        hdy = 10.

        glDisable(GL_LIGHTING)
        glDisable(GL_CULL_FACE)
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0,hx[1],0,hy[1],0,10)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
     
        # Just bin the density
        nbin = 10 
        dhist = np.zeros([nbin])
        binl = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8]
        binu = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
        #p.rho
        # Draw a histogram
        # Compute non-normalised hist
        for i in range(p.n):
            idx = np.where( (p.rho[i] > binl) & (p.rho[i] < binu) )
            dhist[idx] += 1

        for i in range(nbin):
            dhist[i] = (dhist[i] * (binu[i] - binl[i])) / p.n
        glTranslatef(0.0,0.0,-5)
        glColor3f(1.0,0.0,0.5)
        for i in range(nbin):
            glBegin(GL_POLYGON)
            glVertex3f(0.0, 0.0, 0.0)
            glVertex3f(0.0, dhist[i]*hdy, 0.0)
            glVertex3f(hdx, dhist[i]*hdy, 0.0)
            glVertex3f(hdx, 0.0, 0.0)
            glEnd()
            glTranslatef(hdx, 0.0, 0.0)
            #lab = plabel("energy label",loc=(10,10),col=(255,0,255,255) )
            #lab.draw()
  
        glFlush()
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glEnable(GL_LIGHTING)
        glEnable(GL_CULL_FACE)

    def redraw(self,p):
        t = time()
        self.win.dispatch_events()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #self.render_density(p)
        self.hud(p)
        #self.draw_neighbours(p)
        self.draw_particles(p)
        self.density_histogram(p)
        self.draw_box()
        self.timing['Draw time'] = time() - t
