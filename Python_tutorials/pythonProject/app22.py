# nested loops

import os
os.system("clear")

for x in range(4):
    for y in range(3): # for each iteration, the inner loop performs a loop
        print(f"({x}, {y})")