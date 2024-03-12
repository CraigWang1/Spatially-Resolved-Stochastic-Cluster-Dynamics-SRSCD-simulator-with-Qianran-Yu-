VOLUME = 1e-23
LENGTH = 20
SURFACE_LEN = 0.544
SURFACE_VOLUME = VOLUME / LENGTH * SURFACE_LEN

num = 200
num_surface = round(num / VOLUME * SURFACE_VOLUME)


with open('restart.txt', 'w') as f:
    f.write("step = 0\n")
    f.write("time = 0.0\n")
    f.write("fluenceH = 0.0\n")
    f.write(f"object 1    {num_surface}")
    for i in range(100):
        f.write(f"    {num}")

    # Now write hydrogen release
    # f.write("\n")
    # f.write(f"object 1    0")
    # f.write(f"    {1200 * 100}")

    # for i in range(50):
    #     f.write("    0")


    # for i in range(49):
    #     f.write("    0")
