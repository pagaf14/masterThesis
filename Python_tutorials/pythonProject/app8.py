course = "Python for Beginners"
print(len(course))
print(course.upper()) # converts the string into uppercase
print(course) # still in its original form
print(course.lower())
print(course.find("P")) # 0 ==> the index of P in the string is 0
print(course.find("o")) # 4 ==> the index of o in the string is 4
print(course.find("O")) # -1 ==> no character like this in the string
print(course.find("Beginners")) # 11 ==> this sequence of characters starts at index 11
print(course.replace("Beginners", "Absolute Masters"))
# if it doesn't find the old sequence in the original string, it gives back
# the same original string
print("Python" in course) # this expression produces a boolean value, like True or False
