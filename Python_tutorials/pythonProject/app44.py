import os
os.system("clear")

class Person:
    def __init__(self, name):
        self.name = name
        
    def talk(self):
        print(f"Hi, I am {self.name}!")


matteo = Person("Matteo Paganelli")
print(matteo.name)
matteo.talk()