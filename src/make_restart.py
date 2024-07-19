VOLUME = 1e-17  #cm^3
LENGTH = 20     #nm
SURFACE_LEN = 0.544  #nm
SURFACE_VOLUME = VOLUME / LENGTH * SURFACE_LEN + VOLUME  #cm^3
POINTS = 300

num = 3
num_surface = round(num / VOLUME * SURFACE_VOLUME)


with open('restart.txt', 'w') as f:
    f.write("step = 0\n")
    f.write("time = 0.0\n")
    f.write("fluenceH = 0.0\n")
    f.write(f"object -1000000    {num_surface}")
    for i in range(POINTS):
        f.write(f"    {num}")
    f.write("\n")

    """
    f.write(f"object 50000000    1")
    for i in range(100):
        f.write(f"    0")
    f.write("\n")

    f.write(f"object -510000100    1")
    for i in range(100):
        f.write(f"    0")
    """
    

    # Now write hydrogen release
    # f.write("\n")
    # f.write(f"object 1    0")
    # f.write(f"    {1200 * 100}")

    # for i in range(50):
    #     f.write("    0")


    # for i in range(49):
    #     f.write("    0")
