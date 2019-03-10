#!/usr/bin/env python3

# try to replicate a model wave using a spring force simulation



### CONFIG

import numpy as np
import math

SCREEN_SIZE = 800, 400
SAMPLE_RATE = 44100
FPS = 60
ORIGIN_COLOR =  (63, 63, 63) # color of the line through the middle of the screen
MODEL_COLOR =   (255, 0, 0) # color of the perfect model
REPLICA_COLOR = (0, 255, 0) # color of the replica
DOMAIN = np.arange(0, 4, 1/SAMPLE_RATE) # domain of the plot
DISPLAY_TIME_WINDOW = .05 # display this many seconds worth of samples
SPRING_STRENGTH = lambda t: 20000000 * (math.sin(t * 2 * math.pi) * 0.75 + 1)
SPRING_DAMPING = lambda t: 100000000 * (math.cos(t * 2 * math.pi) * 0.75 + 1)
SYNC_SCOPE = True
MODEL_FILE_NAME = 'out_model.wav'
SIMULATION_FILE_NAME = 'out_simulation.wav'



### WAVE STUFF

class Simulation:
    ''' simulate the replica numerically '''
    def __init__(self, model):
        self.model = model

    def integrate(self, ppx, px, a, dt):
        ''' integrate x+1 give x and x-1, acceleration and delta time '''
        return 2 * px - ppx + a * dt ** 2

    def spring(self, pa, a, pb, b, t, k=SPRING_STRENGTH, c=SPRING_DAMPING):
        ''' return the acceleration due to a spring between a and b on a '''
        if callable(k): k = k(t)
        if callable(c): c = c(t)
        f = (b - a) * k # restoring force
        va = a - pa # a velocity
        vb = b - pb # b velocity
        v = vb - va # b relative velocity
        d = c * v # damping force
        return f + d # ignoring mass differences

    def plot(self, domain):
        ''' run the simulation and plot the results '''

        # fill the plot with the first two points in the model
        # this is enough information to get started
        # including initial position and velocity
        plot = self.model[:2]

        # run through the plot domain
        for i in range(2, len(domain)):
            t = domain[i]
            dt = t - domain[i - 1]
            ppx = plot[-2]
            px = plot[-1]
            ppm = self.model[i - 1]
            pm = self.model[i]
            a = self.spring(ppx, px, ppm, pm, t)
            plot.append(self.integrate(ppx, px, a, dt))

        # return the plot
        return plot

def flatten(t):
    ''' flatten a tuple '''
    return sum(t, tuple())

def plot(func, domain):
    ''' plot a function given a domain '''
    return list(map(func, domain))

def resample_plot(plot, size, offset, window_size):
    ''' down or upscale a plot without interpolation '''
    return list(map(lambda i: plot[int(offset + i / size * window_size) % len(plot)], range(size)))

def sync(plot, offset):
    ''' return the position of the first sample that croses the origin '''
    last = None
    check = None
    while not last or last == check:
        last = check
        check = plot[offset % len(plot)] < 0
        offset += 1
    return offset

def save_plot(plot, filename):
    ''' write a plot to a .wav file '''
    with wave.open(filename, 'w') as w:
        w.setnchannels(1)
        w.setsampwidth(1)
        w.setframerate(SAMPLE_RATE)
        for sample in plot:
            w.writeframesraw(bytes([max(0, min(255, int(sample * 255)))]))

def load_plot(filename):
    ''' read a plot from a .wav file '''
    def read_frame(frame):
        ''' convert a frame into a plot sample '''
        frame = w.readframes(1)
        # get the samples for each channel
        def get_sample(i):
            width = w.getsampwidth()
            start = i * width
            end = start + width
            sample = int.from_bytes(frame[start:end], byteorder='little', signed=True)
            return sample / 2 ** (8 * width)
        samples = list(map(get_sample, range(w.getnchannels())))
        # and average them for mononess bro
        sample = sum(samples) / len(samples)
        return sample
    with wave.open(filename, 'r') as w:
        return list(map(read_frame, range(w.getnframes())))



### PLOTTING

import sys
import wave

# see if a wav file was specified for loading!
if len(sys.argv) > 1:
    fn = sys.argv[1]
    print(f'loading model from {fn}...')
    model_plot = load_plot(fn)
else:
    print('plottng model...')
    model_func = lambda t: (1 - t * 110 % 1 * 2) * .4 # sawtooth waveform
    model_plot = plot(model_func, DOMAIN)

print('plotting simulation...')
replica_plot = Simulation(model_plot).plot(DOMAIN)

print(f'writing model result to {MODEL_FILE_NAME}...')
save_plot(model_plot, MODEL_FILE_NAME)

print(f'writing simulation result to {SIMULATION_FILE_NAME}...')
save_plot(replica_plot, SIMULATION_FILE_NAME)



### PYGLET STUFF

import pyglet
import time
window = pyglet.window.Window(*SCREEN_SIZE)

@window.event
def on_draw():
    ''' draw the screen '''

    def vertex(i, plot):
        ''' return a vertex for plotting a point from a plot on screen '''
        h = SCREEN_SIZE[1] / 2
        return i / len(plot) * SCREEN_SIZE[0], plot[i] * h + h

    def draw_plot(plot, color):
        ''' draw a plot to the screen '''
        vertices = flatten(tuple(map(lambda i: vertex(i, plot), range(len(plot)))))
        pyglet.graphics.draw(len(plot), pyglet.gl.GL_LINE_STRIP,
            ('v2f', vertices),
            ('c3B', color * len(plot)),
            )

    # clear the screen first
    window.clear()

    # draw the origin line
    w, y = SCREEN_SIZE[0], SCREEN_SIZE[1] / 2
    pyglet.graphics.draw(2, pyglet.gl.GL_LINE_STRIP,
        ('v2f', (0, y, w, y)),
        ('c3B', ORIGIN_COLOR * 2),
        )

    # calculate the offset
    t = time.time()
    offset = int(t * SAMPLE_RATE)
    window_size = DISPLAY_TIME_WINDOW * SAMPLE_RATE

    # optionally syncronize the scope to the origin of the waveform
    if SYNC_SCOPE: offset = sync(model_plot, offset)

    # draw the model plot
    model_plot_render = resample_plot(model_plot, SCREEN_SIZE[1], offset, window_size)
    draw_plot(model_plot_render, MODEL_COLOR)

    # draw the replica plot
    replica_plot_render = resample_plot(replica_plot, SCREEN_SIZE[1], offset, window_size)
    draw_plot(replica_plot_render, REPLICA_COLOR)

def update(dt):
    pass

print('running...')
pyglet.clock.schedule_interval(update, 1 / FPS)
pyglet.app.run()
