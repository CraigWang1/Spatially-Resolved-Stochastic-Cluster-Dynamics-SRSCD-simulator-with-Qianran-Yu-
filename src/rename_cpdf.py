import os

for i in range(100, 0, -1):
    os.rename(f"cpdf{i}.txt", f"cpdf{i+1}.txt")
