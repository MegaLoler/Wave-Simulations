#!/usr/bin/env python3

# simulate a transverse wave on a string in 2d space
# using spring forces
# this will be a live numerical solution



### CONFIG

SCREEN_SIZE = 500, 200
STRING_RES = SCREEN_SIZE[0]
FPS = 60
SWEEPS = 20 # screen widths per second
ORIGIN_COLOR = (63, 63, 63) # color of the line through the middle of the screen
STRING_COLOR = (255, 0, 0) # color of the string
SPRING_STRENGTH = .5
SPRING_DAMPING = .05



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

    def spring(self, a, pa, b, pb, k=SPRING_STRENGTH, c=SPRING_DAMPING):
        ''' return the acceleration due to a spring between a and b on a '''
        f = (b - a) * k * .5 # restoring force
        va = a - pa # a velocity
        vb = b - pb # b velocity
        v = vb - va # b relative velocity
        d = c * v # damping force
        return f + d # ignoring mass differences

    def update(self, t, dt):
        ''' update the simulation '''

        # backup the current state
        current_state = list(self.points)

        # update the rest of the displacement points
        new_state = [self.function(t)]
        for i, p, x in zip(range(len(self.points)), self.points_p, self.points):
            v = x - p # velocity
            a = 0 # acceleration
            if i < len(self.points) - 1:
                r = self.points[i + 1]
                rp = self.points_p[i + 1]
                a += self.spring(x, p, r, rp)
            if i > 0:
                l = self.points[i - 1]
                lp = self.points_p[i - 1]
                a += self.spring(x, p, l, lp)
            if i > 0 and i < len(self.points) - 1:
                new_state.append(x + v + a)# * dt ** 2)
        new_state.append(0)

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
    for i in range(SWEEPS):
        simulation.update(time.time() - start_time, dt / SWEEPS)

start_time = time.time()

#wave_func = lambda t: math.sin(t * 2 * math.pi) * 0.5
#wave_func = lambda t: 0.5 if t > 1 and t < 1.5 else 0
wave_func = lambda t: ((t % 1) * 2 - 1) * 0.25
simulation = String(STRING_RES, wave_func)

pyglet.clock.schedule_interval(update, 1 / FPS)
pyglet.app.run()
