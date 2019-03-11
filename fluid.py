#!/usr/bin/env python3

# 2d fluid simulation using lattice boltzmann algorithm



### CONFIG

SCREEN_SIZE = 400, 400
FPS = 60 # desired number of animation frames to render per second
SWEEPS = 1 # number of simulation steps per frame
DIMENSIONS = 50, 50 # dimensions of the simulation lattice
VISCOSITY = 0.02
MOUSE_SENSITIVITY = 0.05
DEBUG = False



### CONSTANTS

OMEGA = 1 / (3 * VISCOSITY + .5) # reciprocal relaxation time
v_4_9 = 4 / 9



### IMPORTS

import itertools
import numpy as np



### CLASSES

class Node:
    ''' represents a lattice node '''

    def __init__(self, rho=1, ux=0, uy=0):
        self.invalidate_cache()
        self.set_equilibrium(ux, uy, rho)

    def invalidate_cache(self):
        ''' invalidate any cached data, meaning it has to be recalculated '''
        self.cache_rho = None

    @property
    def densities(self):
        ''' return a sequence of the discretized densities '''
        return self.c, self.n, self.s, self.e, self.w, self.nw, self.ne, self.sw, self.se

    @property
    def ux(self):
        ''' return the x component of this node's velocity '''
        return (self.e + self.ne + self.se - self.w - self.nw - self.sw) / self.rho

    @property
    def uy(self):
        ''' return the y component of this node's velocity '''
        return (self.n + self.nw + self.ne - self.s - self.sw - self.se) / self.rho

    @property
    def u(self):
        ''' return this node's velocity '''
        return self.ux, self.uy

    @u.setter
    def u(self, value):
        ''' set the macroscopic velocity of this node '''
        ux, uy = value
        self.set_equilibrium(ux, uy, self.rho)

    @property
    def rho(self):
        ''' return the macroscopic density of this node '''
        if not self.cache_rho:
            self.cache_rho = sum(self.densities)
        return self.cache_rho

    @rho.setter
    def rho(self, value):
        ''' set the macroscopic density of this node '''
        self.set_equilibrium(self.ux, self.uy, value)

    def set_equilibrium(self, ux, uy, rho):
        ''' set the macroscopic velocity and density according to the equilibrium distribution '''
        # TODO: actually understand this stuff lol
        # rn its just taken from https://physics.weber.edu/schroeder/fluids/
        v_1_9_rho = rho / 9
        v_1_36_rho = rho / 36
        ux3 = 3 * ux
        uy3 = 3 * uy
        ux2 = ux ** 2
        uy2 = uy ** 2
        uxuy2 = 2 * ux * uy
        u2 = ux2 + uy2
        u215 = 1.5 * u2
        self.c  = v_4_9 * rho * (1                                  - u215)
        self.e  = v_1_9_rho   * (1 + ux3       + 4.5 * ux2          - u215)
        self.w  = v_1_9_rho   * (1 - ux3       + 4.5 * ux2          - u215)
        self.n  = v_1_9_rho   * (1 + uy3       + 4.5 * uy2          - u215)
        self.s  = v_1_9_rho   * (1 - uy3       + 4.5 * uy2          - u215)
        self.ne = v_1_36_rho  * (1 + ux3 + uy3 + 4.5 * (u2 + uxuy2) - u215)
        self.se = v_1_36_rho  * (1 + ux3 - uy3 + 4.5 * (u2 - uxuy2) - u215)
        self.nw = v_1_36_rho  * (1 - ux3 + uy3 + 4.5 * (u2 - uxuy2) - u215)
        self.sw = v_1_36_rho  * (1 - ux3 - uy3 + 4.5 * (u2 + uxuy2) - u215)
        self.cache_rho = rho

    def collide(self):
        ''' update the distribution of mass in this node '''
        # TODO: actually understand this stuff lol
        # rn its just taken from https://physics.weber.edu/schroeder/fluids/
        rho = self.rho
        ux = self.ux
        uy = self.uy
        v_1_9_rho = rho / 9
        v_1_36_rho = rho / 36
        ux3 = 3 * ux
        uy3 = 3 * uy
        ux2 = ux ** 2
        uy2 = uy ** 2
        uxuy2 = 2 * ux * uy
        u2 = ux2 + uy2
        u215 = 1.5 * u2
        ux2_4_5 = 4.5 * ux2
        uy2_4_5 = 4.5 * uy2
        u2_p_uxuy2_4_5 = 4.5 * (u2 + uxuy2)
        u2_m_uxuy2_4_5 = 4.5 * (u2 - uxuy2)
        self.c  += OMEGA * (v_4_9 * rho * (1                              - u215) - self.c)
        self.e  += OMEGA * (v_1_9_rho   * (1 + ux3       + ux2_4_5        - u215) - self.e)
        self.w  += OMEGA * (v_1_9_rho   * (1 - ux3       + ux2_4_5        - u215) - self.w)
        self.n  += OMEGA * (v_1_9_rho   * (1 + uy3       + uy2_4_5        - u215) - self.n)
        self.s  += OMEGA * (v_1_9_rho   * (1 - uy3       + uy2_4_5        - u215) - self.s)
        self.ne += OMEGA * (v_1_36_rho  * (1 + ux3 + uy3 + u2_p_uxuy2_4_5 - u215) - self.ne)
        self.se += OMEGA * (v_1_36_rho  * (1 + ux3 - uy3 + u2_m_uxuy2_4_5 - u215) - self.se)
        self.nw += OMEGA * (v_1_36_rho  * (1 - ux3 + uy3 + u2_m_uxuy2_4_5 - u215) - self.nw)
        self.sw += OMEGA * (v_1_36_rho  * (1 - ux3 - uy3 + u2_p_uxuy2_4_5 - u215) - self.sw)
        self.invalidate_cache()

class Lattice:
    ''' represents a lattice of nodes '''

    def __init__(self, dim, fill=Node):
        ''' recursively generate an n-dimensional lattice filled with fill '''

        def make_fill():
            ''' allow fill to be a function or a value '''
            return fill() if callable(fill) else fill

        self.content = list(map(lambda i: Lattice(dim[1:], fill), range(dim[0]))) if dim else make_fill()

    def __len__(self):
        ''' return the number of items in this lattice '''
        return np.prod(self.dimensions)

    def __getitem__(self, coords):
        ''' return the item in the lattice given its coordiates '''
        return self.content[coords[0] % self.dimensions[0]][coords[1:]] if coords else self.content

    def __iter__(self):
        ''' return an iterator to iterate linearly through all nodes '''
        return itertools.chain(*self.content) if self.dimensionality else iter([self.content])

    @property
    def count(self):
        ''' return the length of the top level dimension '''
        return len(self.content) if isinstance(self.content, list) else 0

    @property
    def dimensions(self):
        ''' return the dimensions of this lattice '''
        return (self.count,) + self.content[0].dimensions if self.count else ()

    @property
    def dimensionality(self):
        ''' return the dimensionality of this lattice (2d, 3d, 4d, etc) '''
        return len(self.dimensions)

    def sum(self, func):
        ''' return the sum of applying a function to each node '''
        return sum(np.array(list(map(func, self))))

    def average(self, func):
        ''' return the average of applying a function to each node '''
        return self.sum(func) / len(self)

class Simulation:
    ''' represents a fluid simulation state '''

    def __init__(self, dimensions=DIMENSIONS):
        # "double buffering"
        self.lattice = Lattice(dimensions)
        self.buffer = Lattice(dimensions)

        # for speed reasons cache the list of coordinates and nodes for iteration later
        self.coords = tuple(np.ndindex(dimensions))
        self.lattice_nodes = tuple(self.lattice)
        self.buffer_nodes = tuple(self.buffer)

        # cache the neighbor cells as well
        self.lattice_neighbors = self.cache_neighbors(self.lattice)
        self.buffer_neighbors = self.cache_neighbors(self.buffer)

    def cache_neighbors(self, lattice):
        ''' return a tuple of objects containing neighboring cells '''

        class Neighborhood:
            ''' represent all the neighbors for a node '''
            def __init__(self, coords):
                x, y = coords
                self.c = lattice[x, y]
                self.n = lattice[x, y - 1]
                self.s = lattice[x, y + 1]
                self.e = lattice[x + 1, y]
                self.w = lattice[x - 1, y]
                self.nw = lattice[x - 1, y - 1]
                self.ne = lattice[x + 1, y - 1]
                self.sw = lattice[x - 1, y + 1]
                self.se = lattice[x + 1, y + 1]

        return tuple(map(Neighborhood, self.coords))

    def step(self):
        ''' perform a single step of the simulation '''
        self.collide()
        self.stream()

    def collide(self):
        ''' perform the inner-node collisions '''
        for node in self.lattice_nodes:
            node.collide()

    def stream(self):
        ''' move the mass between nodes according to their velocities '''

        # loop through each node of the new lattice, calculating the values from the old one
        for node, neighborhood, coords in zip(self.buffer_nodes, self.lattice_neighbors, self.coords):
            # TODO: make this..... n dimensional
            # FOR NOW... assuming 2d

            # get the neighbors
            c = neighborhood.c
            n = neighborhood.n
            s = neighborhood.s
            e = neighborhood.e
            w = neighborhood.w
            nw = neighborhood.nw
            ne = neighborhood.ne
            sw = neighborhood.sw
            se = neighborhood.se

            # move the densities
            node.c = c.c
            node.n = s.n
            node.s = n.s
            node.e = w.e
            node.w = e.w
            node.nw = se.nw
            node.ne = sw.ne
            node.sw = ne.sw
            node.se = nw.se

            # invalidate the cache of that node
            node.invalidate_cache()

        # set the current lattice to the new one
        # which amounts to swapping 'buffers'
        self.old_lattice = self.lattice
        self.old_nodes = self.lattice_nodes
        self.old_neighbors = self.lattice_neighbors
        self.lattice = self.buffer
        self.lattice_nodes = self.buffer_nodes
        self.lattice_neighbors = self.buffer_neighbors
        self.buffer = self.old_lattice
        self.buffer_nodes = self.old_nodes
        self.buffer_neighbors = self.old_neighbors

    @property
    def mass(self):
        ''' return the total mass in the system '''
        return self.lattice.sum(lambda n: n.rho)

    @property
    def velocity(self):
        ''' return the total mass in the system '''
        return self.lattice.average(lambda n: n.u)

    def draw(self, buffer_size):
        ''' render the state of the simulation, returning raw image bytes data '''

        # allocate the output buffer
        image_buffer = bytearray(buffer_size)

        # loop through lattice coordinates and collect vertices and colors
        i = 0
        for node in self.lattice_nodes:
            # append the color for it
            # TODO: better, more general visualization
            value = max(0, min(255, int(node.rho * 255) - 200))
            r = g = b = value

            # add the color to the image buffer
            image_buffer[i] = r
            image_buffer[i + 1] = g
            image_buffer[i + 2] = b
            image_buffer[i + 3] = 255
            i += 4

        # return the buffer
        return bytes(image_buffer)



### MAIN

def main():
    ''' run a simulation '''

    # create the simulation
    simulation = Simulation()

    # cause an initial disturbance in the middle of the screen
    #for x, y in np.ndindex(3, 3):
    #    simulation.lattice[x + DIMENSIONS[0] // 2, y + DIMENSIONS[1] // 2].rho = 5

    # use pyglet to display the simulation in real time
    import pyglet
    import pyglet.gl
    import pyglet.font
    import pyglet.window

    # set opengl parameters
    pyglet.gl.glEnable(pyglet.gl.GL_TEXTURE_2D)
    pyglet.gl.glTexParameteri(pyglet.gl.GL_TEXTURE_2D, pyglet.gl.GL_TEXTURE_MAG_FILTER, pyglet.gl.GL_NEAREST)
    
    # create pyglet window
    window = pyglet.window.Window(*SCREEN_SIZE)

    # create the drawing surface
    surface = pyglet.image.SolidColorImagePattern((0, 0, 0, 0,),).create_image(*DIMENSIONS)
    surface_size = np.prod(DIMENSIONS) * 4

    # keep track of actual fps
    fps_display = pyglet.clock.ClockDisplay(font=pyglet.font.load('Mono', 8, bold=True), color=(1,1,0,.5))

    @window.event
    def on_draw():
        ''' draw the screen '''

        # clear the screen first
        window.clear()

        # render the simulation to the drawing surface
        surface.set_data('RGBA', DIMENSIONS[0] * 4, simulation.draw(surface_size))

        # display the rendered image on screen
        surface.texture.width, surface.texture.height = SCREEN_SIZE
        surface.blit(0, 0)

        # display fps
        fps_display.draw()

    @window.event
    def on_mouse_drag(x, y, dx, dy, button, modifiers):
        ''' make the mouse able to drag particles '''
        if button == pyglet.window.mouse.LEFT:
            # get the simulation coordinates
            x, y = x / SCREEN_SIZE[0] * DIMENSIONS[0], y / SCREEN_SIZE[1] * DIMENSIONS[1]
            x, y = int(x), int(y)
            node = simulation.lattice[y, x]
            ux, uy = node.u
            dx, dy = dx * MOUSE_SENSITIVITY, dy * MOUSE_SENSITIVITY
            node.u = dy + ux, -dx + uy

    def update(dt):
        ''' update the simulation for the next animation frame '''

        # run a number of simulation steps for this frame
        for i in range(SWEEPS):
            simulation.step()

        # show some useful information
        if DEBUG:
            print(f'system mass     = {simulation.mass}')
            print(f'system velocity = {simulation.velocity}')

    # setup the pyglet update interval
    pyglet.clock.schedule_interval(update, 1 / FPS)

    # run the interface
    pyglet.app.run()

if __name__ == '__main__': main()
