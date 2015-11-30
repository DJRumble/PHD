file = open("datafile.txt", "r")
filedata = file.read()
file.close()

newLine = "Your new line of data with the time stamp goes here.\n" + 
filedata
file = open("datafile.txt", "w")
file.write(newLine)
file.close()
