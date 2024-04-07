import os
os.system('clear')

message = input("> ")
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
    output += emojis.get(word, word) + " "  # if the word is not a key in the dictionary use the original word
print(output)