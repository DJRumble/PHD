file = open("20130306_extraneous_cores_aquila.txt", "r")
filedata = file.read()
file.close()

newLine = "Your new line of data with the time stamp goes here.\n" + filedata
file = open("20130306_extraneous_cores_aquila.txt", "w")
file.write(newLine)
file.close()
