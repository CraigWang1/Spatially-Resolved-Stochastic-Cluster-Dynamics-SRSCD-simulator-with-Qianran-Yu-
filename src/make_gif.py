import os
import imageio

images = []
for f in sorted(os.listdir('images')):
	images.append(imageio.imread('images/' + f))
imageio.mimsave('saturated.gif', images, loop=0, duration=100)