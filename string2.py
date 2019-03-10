#!/usr/bin/env python3

# simulate a transverse wave on a string in 2d space
# this will be a live numerical solution



### CONFIG

SCREEN_SIZE = 500, 200
STRING_RES = SCREEN_SIZE[0]
WAVE_VELOCITY = 0.25 # in units of the length of the string (width of screen) per second
FPS = 60
ORIGIN_COLOR = (63, 63, 63) # color of the line through the middle of the screen
STRING_COLOR = (255, 0, 0) # color of the string



### STRING STUFF

def flatten(t):
    ''' flatten a tuple '''
    return sum(t, tuple())

class String:
    ''' represents a string to simulate '''

    def __init__(self, resolution, function):
        # these are the displacement points
        self.points = [0] * resolution

        # and this is the previous state of them
        self.points_p = list(self.points)

        # the function to fixate to the end of the string
        self.function = function

    def draw(self, window, color=STRING_COLOR):
        ''' draw the string onto a pyglet window '''
        def vertex(i):
            ''' make a vertex for a displacement point '''
            h = SCREEN_SIZE[1] / 2
            return i / len(self.points) * SCREEN_SIZE[0], self.points[i] * h + h
        res = len(self.points)
        vertices = list(map(vertex, range(res)))
        pyglet.graphics.draw(res, pyglet.gl.GL_LINE_STRIP,
            ('v2f', flatten(vertices)),
            ('c3B', color * res),
            )

    def spring(self, a, ap, b, bp):
        ''' return acceleration due to spring '''
        

    def update(self, t, dt):
        ''' update the simulation '''

        # backup the current state
        current_state = list(self.points)

        # set the end to the wave function
        self.points[0] = self.function(t)

        # update the rest of the displacement points
        new_state = list()
        for i, p, x in zip(len(self.points), self.points_p, self.points):
            v = x - p # velocity
            a = 0 # acceleration
            if i < len(self.points) - 1:
                r = self.points[i + 1]
                rp = self.points_p[i + 1]
                a += self.spring(x, p, r, rp, dt)
            if i > 0:
                l = self.points[i - 1]
                lp = self.points_p[i - 1]
                a += self.spring(x, p, l, lp, dt)
            new_state.append(x + v + a * dt ** 2)

        # set the new state
        self.points = new_state
        self.points_p = current_state



### PYGLET STUFF

import time
import math
import pyglet
window = pyglet.window.Window(*SCREEN_SIZE)

@window.event
def on_draw():
    ''' draw the screen '''

    # clear the screen first
    window.clear()

    # draw the origin line
    w, y = SCREEN_SIZE[0], SCREEN_SIZE[1] / 2
    pyglet.graphics.draw(2, pyglet.gl.GL_LINE_STRIP,
        ('v2f', (0, y, w, y)),
        ('c3B', ORIGIN_COLOR * 2),
        )

    # draw the string
    simulation.draw(window)

def update(dt):
    simulation.update(time.time(), dt)

wave_func = lambda t: math.sin(t * 2 * math.pi)
simulation = String(STRING_RES, wave_func)

pyglet.clock.schedule_interval(update, 1 / FPS)
pyglet.app.run()
