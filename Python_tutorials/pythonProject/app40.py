# re-organize app36.py in a function

import os
os.system('clear')

def emoji_converter(message):
    words = message.split(' ')
    # print(words)
    emojis = {
        ":)": "😁",
        ":D": "😃",
        ";)": "😉",
        ":(": "😔"
    }
    output = ""
    for word in words:
        output += emojis.get(word, word) + " "
    return output


message = input("> ") # not in the function because input can come from different sources (graph-code user interface)
print(emoji_converter(message)) # out of function cause what to do with output is different from case to case