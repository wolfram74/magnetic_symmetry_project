import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import math
from PIL import Image

class collage_builder():

    def __init__(self):
        self.num_graphs=0
        self.files_list = []

    def add_solution(self, positions, rotations, u, tau, nature):
        self.num_graphs += 1
        n = len(positions)
        circles = [plt.Circle((positions[i][0], positions[i][1]), 0.5, color=vec_colour(rotations[i]), ec = 'k') for i in range(len(positions))]

        fig, ax = plt.subplots(figsize=(5, 5))

        ax.set_aspect('equal')
        ax.set_ylim(-2.5, 2.5)
        ax.set_xlim(-2.5, 2.5)

        for c in circles:
            ax.add_artist(c)

        for i in range(n):
            ax.quiver([positions[i][0]], [positions[i][1]], [rotations[i][0]], [rotations[i][1]], units='xy', scale=2, color=vec_colour2(rotations[i]), zorder=10, width=0.1)
        
        dipole_moment = np.linalg.norm(sum(rotations))
        ax.text(-2.4, 2.1, 'sol. no. ' + str(self.num_graphs) + ', nature = ' + nature)
        ax.text(-2.4, -2.1, 'U='+str('%.6g'%u)+'  m='+str('%.6g'%dipole_moment)+'  tau^2='+str('%.6g'%tau), color='k')
        plt.savefig('temp'+str(self.num_graphs)+'.jpg')
        self.files_list.append('temp'+str(self.num_graphs)+'.jpg')


    def create_collage(self, filename, width=2000, height=3000):
        cols = 6
        rows = 9
        thumbnail_width = width//cols
        thumbnail_height = height//rows
        size = thumbnail_width, thumbnail_height
        new_im = Image.new('RGB', (width, height))
        ims = []
        for p in self.files_list:
            im = Image.open(p)
            im.thumbnail(size)
            ims.append(im)
        i = 0
        x = 0
        y = 0
        for col in range(cols):
            for row in range(rows):
                if i == len(self.files_list):
                    break
                print(i, x, y)
                new_im.paste(ims[i], (x, y))
                i += 1
                y += thumbnail_height
            x += thumbnail_width
            y = 0

        new_im.save(filename)


def vec_colour(v):
    theta = math.atan2(v[1], v[0])
    return ((1+math.cos(theta))/2,0.0,(1+math.sin(theta))/2)

def vec_colour2(v):
    theta = math.atan2(v[1], v[0])
    return ((1-math.cos(theta))/2,1.0,(1-math.sin(theta))/2)
