import os
os.system("clear")

command = ""
started = False

while True: # != means different; equivalent: command.upper ... QUIT
    command = input("> ").lower() # otherwise you would have to write every time command.lower()
    if command == "start":
        if started:
            print("Car is already started!")
        else:
            started = True
            print("Car started ... Ready to go!")
    elif command == "stop":
        if not started:
            print("Car is already st.opped!")
        else:
            started = False
            print("Car stopped")
    elif command == "help":
        print("""
start - to start the car
stop - to stop the car
quit - to quit
          """)
    elif command == "quit":
        break
    else:
        print("Sorry, I don't understand that")