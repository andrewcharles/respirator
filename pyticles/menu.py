import time

import pyglet
from pyglet import graphics
from pyglet import text
from pyglet import event
from pyglet.gl import *

class SimpleMenu(event.EventDispatcher, graphics.Batch):
    def __init__(self, parent, x, y, options):
        super(SimpleMenu, self).__init__()

        self.parent = parent
        self.x, self.y = x, y

        self.options = []
        y = y
        self.width = 0
        for option in options:
            l = text.Label(option, x=x, y=y, color=(0, 0, 0, 255), batch=self)
            self.options.append(l)
            self.width = max(self.width, l.content_width)
            y -= l.content_height

        self.height = abs(self.y - y)

        # adjust menu position to make sure whole menu is visible
        if self.x < 0:
            self.x = 0
        if self.x + self.width> parent.width:
            self.x = parent.width - self.width
        if self.y - self.height < 0:
            self.y = self.height
        if self.y > parent.height:
            self.y = parent.height

        y = self.y
        for option in self.options:
            option.x = self.x
            option.y = y
            y -= l.content_height

        self.created_time = time.time()

        parent.push_handlers(self)

    def on_mouse_motion(self, x, y, dx, dy):
        ret = event.EVENT_UNHANDLED
        for option in self.options:
            if (option.x < x < option.x + self.width and
                    option.y - option.content_height < y < option.y):
                option.color = (200, 0, 0, 255)
                ret = event.EVENT_HANDLED
            elif option.color != (0, 0, 0, 255):
                option.color = (0, 0, 0, 255)
        return ret

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        return self.on_mouse_motion(x, y, dx, dy)

    def on_mouse_press(self, x, y, buttons, modifiers):
        # determine action to take
        selected_option = None
        for option in self.options:
            if (option.x < x < option.x + self.width and
                    option.y - option.content_height < y < option.y):
                selected_option = option.text
                break

        self.parent.remove_handlers(self)
        self.options = None
        self.dispatch_event('on_close', selected_option)
        return event.EVENT_HANDLED

    def on_mouse_release(self, x, y, button, modifiers):
        if time.time() - self.created_time > .5:
            return self.on_mouse_press(x, y, button, modifiers)
        return event.EVENT_UNHANDLED

    def draw(self):
        glColor4f(255, 255, 255, 255)
        glRectf(self.x, self.y, self.x + self.width, self.y - self.height)
        super(SimpleMenu, self).draw()

SimpleMenu.register_event_type('on_close')

if __name__ == '__main__':
    win = pyglet.window.Window(height=400)
    m = None

    @win.event
    def on_draw():
        win.clear()
        if m: m.draw()

    @win.event
    def on_mouse_press(x, y, button, modifiers):
        global m
        m = SimpleMenu(win, x, y, ['New', 'Load', 'Save', 'Quit'])
        @m.event
        def on_close(selected_option):
            print selected_option
            global m
            m = None
        return event.EVENT_HANDLED

    pyglet.app.run()

