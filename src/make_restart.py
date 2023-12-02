VOLUME = 1e-23
LENGTH = 20
SURFACE_LEN = 0.544
SURFACE_VOLUME = VOLUME / LENGTH * SURFACE_LEN

num = 100
num_surface = round(num / VOLUME * SURFACE_VOLUME)


with open('restart.txt', 'w') as f:
    f.write("step = 0\n")
    f.write("time = 0.0\n")
    f.write("fluenceH = 0.0\n")
    f.write(f"object -1000000    {num_surface}")
    for i in range(100):
        f.write(f"    {num}")
