#!/usr/bin/env python3

# simulate a transverse wave on a string in 2d space
# this will be an analytical solution
# the state of the string will be computed as a function of time



### CONFIG

SCREEN_SIZE = 500, 200
STRING_RES = SCREEN_SIZE[0]
WAVE_VELOCITY = 0.25 # in units of the length of the string (width of screen) per second
FPS = 60
ORIGIN_COLOR =    (63, 63, 63) # color of the line through the middle of the screen
COMPONENT_COLOR = (0, 63, 127) # color of individual strings
SUM_COLOR =       (255, 0, 0) # color of sum of strings



### WAVE STUFF

import time
import math

def wave_func(function, amp=1, freq=1, phase=0, fperiod=1, color=COMPONENT_COLOR):
    ''' return a function that returns amplitude given position and time given a general function

    amplitude is the multiplication factor of the result
    frequency scales the input to the given function
    phase offsets the input to the given function
    fperiod is the period of the given function, used to normalize the function's period
    '''

    def func(x, t):
        # see if any of the parameters are functions!!!
        nonlocal amp, freq, phase
        a = amp(t)   if callable(amp)   else amp
        f = freq(t)  if callable(freq)  else freq
        p = phase(t) if callable(phase) else phase
        return a * function((t + p + x / WAVE_VELOCITY) * f * fperiod)

    return func, color

def wave_sum(*funcs):
    ''' return a function that is the sum of multiple wave functions '''
    return lambda x, t: sum(map(lambda f: f(x, t), funcs))

def flatten(t):
    ''' flatten a tuple '''
    return sum(t, tuple())

def draw_string(window, wave_function, time, color=(255, 255, 255)):
    ''' draw the string in its current state '''

    def point(i):
        ''' return the on-screen coordinates for a point on the string '''
        displacement = wave_function(i / STRING_RES, time)# / len(wave_functions)
        h = SCREEN_SIZE[1] / 2
        x, y = i * SCREEN_SIZE[0] / STRING_RES, h + displacement * h
        return x, y

    length = STRING_RES + 1
    vertices = tuple(map(point, range(length)))
    pyglet.graphics.draw(length, pyglet.gl.GL_LINE_STRIP,
        ('v2f', flatten(vertices)),
        ('c3B', color * length),
        )

def set_wave_functions(*funcs):
    ''' set the wave functions to be drawn as strings '''
    global wave_functions, wave_function
    wave_functions = funcs
    wave_function = wave_sum(*tuple(map(lambda f: f[0], funcs)))

# define the functions to draw here dude
set_wave_functions(
        *map(lambda i: wave_func(math.sin, amp=.5/i, freq=i/2, phase=.0, fperiod=2*math.pi),
            range(1, 20))
        )



### PYGLET STUFF

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

    # draw a string for each individual function
    t = time.time()
    for func, color in wave_functions:
        draw_string(window, func, t, color)

    # draw a string for the sum of the functions
    draw_string(window, wave_function, t, SUM_COLOR)

def update(dt):
    pass

pyglet.clock.schedule_interval(update, 1 / FPS)
pyglet.app.run()
