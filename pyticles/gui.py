""" A very minimal library of graphical user interface elements using simple
pyglet primitives.


References:
http://www.python-forum.org/pythonforum/viewtopic.php?f=2&t=11160
Re: Using pyglet to view images
Postby vegaseat on Sat Feb 07, 2009 11:22 am 


"""

import pyglet
from pyglet.window import mouse
from pyglet.gl import *

def nothing():
    pass

class Button:
    """ A button.

        __init__() -- create a button with default values
        color -- button is this color
        tcolor -- text is this color
        mcolor -- button goes this color when mousedover?

    """
    
    def __init__(self
                ,loc = (30,30)
                ,size = (32,32)
                ,color = (0.5,0.5,0.5)
                ,mcolor = (0.9,0.9,0.9)
                ,activate = nothing()
                ,image = None
                ,labeltext = "button"
                ,group = None
                ):
        self.color = color
        self.mcolor = mcolor
        self.x = loc[0]
        self.y = loc[1]
        self.width = size[0] 
        self.height = size[1]
        self.activate = activate
        self.img = image
        self.group = group
        self.label = pyglet.text.Label(labeltext,
            font_name="Arial", 
            font_size=12,color=(220,220,220,244) ,
            x=self.x+5,y=self.y+5,group=group)


    def activate(self):
            print "Hit",self.label.text

    def hit(self,x,y):
        """ Checks if the given coords are in the
            button's hitbox.
        """
        if( (x > self.x) and (y > self.y)
        and (x < self.x + self.width) and (y < self.y + self.width)):
            #self.activate()
            return True
    
    def draw(self):
        """ Draw the button. We push the opengl state so that
            we can draw to the screen plane, and then pop it
            back.
        """

        glClearColor(self.color[0],self.color[1],self.color[2],1.0)
        self.label.draw()

        if self.img:
            self.img.blit(self.x,self.y)
        else:
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
            glBegin(GL_POLYGON)
            glColor3f(self.color[0],self.color[1],self.color[2])
            glVertex2f(self.x,self.y)
            glVertex2f(self.x,self.y+self.height)
            glVertex2f(self.x+self.width,self.y+self.height)
            glVertex2f(self.x+self.width,self.y)
            glVertex2f(self.x,self.y)
            glEnd()
        



