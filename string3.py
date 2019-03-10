#!/usr/bin/env python3

# simulate a transverse wave on a string in 2d space
# this will be a simple progression of energy left and right
# this will be a live numerical solution



### CONFIG

SCREEN_SIZE = 500, 200
STRING_RES = SCREEN_SIZE[0]
FPS = 60
SWEEPS = 1
WAVE_VELOCITY = FPS * SWEEPS / SCREEN_SIZE[0]
ORIGIN_COLOR = (63, 63, 63) # color of the line through the middle of the screen
STRING_COLOR = (255, 0, 0) # color of the string
COLORS = ((255, 0, 0), (0, 255, 0), (0, 0, 255))
REFLECTION = .5



### STRING STUFF

def flatten(t):
    ''' flatten a tuple '''
    return sum(t, tuple())

class String:
    ''' represents a string to simulate '''

    def __init__(self, resolution, function):
        # keep track of energy moving left and right both
        self.left = [0] * resolution
        self.right = [0] * resolution

        # the function to fixate to the end of the string
        self.function = function

    def __len__(self):
        ''' return the resolution of this string '''
        return len(self.left)

    def __getitem__(self, i):
        ''' return the sum energy at a point on this string '''
        return self.left[i] + self.right[i]

    def draw(self, window, color=STRING_COLOR):
        ''' draw the string onto a pyglet window '''
        def vertex(i):
            ''' make a vertex for a displacement point '''
            h = SCREEN_SIZE[1] / 2
            return i / len(self) * SCREEN_SIZE[0], self[i] * h + h
        res = len(self)
        vertices = list(map(vertex, range(res)))
        pyglet.graphics.draw(res, pyglet.gl.GL_LINE_STRIP,
            ('v2f', flatten(vertices)),
            ('c3B', color * res),
            )

    def update(self, t):
        ''' update the simulation '''
        
        # generate the sound source
        source = self.function(t)

        # go round robin
        carry_right = self.right[-1]
        carry_left = self.left[0] + source
        for x in reversed(range(1, len(self))):
            self.right[x] = self.right[x - 1]
        self.right[0] = carry_left * REFLECTION
        for x in range(0, len(self) - 1):
            self.left[x] = self.left[x + 1]
        self.left[len(self) - 1] = carry_right * REFLECTION



### PYGLET STUFF

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
    i = 0
    for simulation in simulations:
        simulation.draw(window, COLORS[i])
        i += 1

def update(dt):
    global time
    for simulation in simulations:
        for i in range(SWEEPS):
            simulation.update(time)
            time += 1 / FPS / SWEEPS

time = 0
simulations = list()
wave_func = lambda t: ((t * WAVE_VELOCITY % 1) * 2 - 1) * 0.5
simulations.append(String(STRING_RES, wave_func))
wave_func = lambda t: math.sin(t * WAVE_VELOCITY * 2 * math.pi) * 0.5
simulations.append(String(STRING_RES, wave_func))
wave_func = lambda t: 0.5 if t > 1 and t < 1.5 else 0
simulations.append(String(STRING_RES, wave_func))

pyglet.clock.schedule_interval(update, 1 / FPS)
pyglet.app.run()
