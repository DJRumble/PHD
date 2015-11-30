import numpy as np

# Define type conversion from c2d table to numpy arrays
type_dict = {}
type_dict['int'] = np.int
type_dict['real'] = np.float32
type_dict['char'] = np.str

type_rev_dict = {}
type_rev_dict[np.int16] = "int"
type_rev_dict[np.int32] = "int"
type_rev_dict[np.int64] = "int"
type_rev_dict[np.float32] = "real"
type_rev_dict[np.float64] = "real"
type_rev_dict[np.str] = "char"
type_rev_dict[np.string_] = "char"
type_rev_dict[str] = "char"


def c2d_read(filename, atpytbl, definition=3):
    '''
    Read a table from a c2d file

    Required Arguments:

        *filename*: [ string ]
            The c2d file to read the table from

    Optional Keyword Arguments:

        *definition*: [ 1 | 2 | 3 ]

            The definition to use to read c2d tables:

            1: any character below a pipe symbol belongs to the
               column on the left, and any characters below the
               first pipe symbol belong to the first column.
            2: any character below a pipe symbol belongs to the
               column on the right.
            3: no characters should be present below the pipe
               symbols (default).
    '''

    if not definition in [1, 2, 3]:
        raise Exception("definition should be one of 1/2/3")

    #reset()

    # Open file for reading
    f = file(filename, 'rb')

    line = f.readline()

    # Read in comments and keywords
    while True:

        char1 = line[0:2]
        char2 = line[1:2]

        if char1 <> '# ':
            break

        if char2==' ' or not '=' in line: # comment
            #add_comment(line[1:])
            pass
        else:          # keyword
            pos = line.index('=')
            key, value = line[1:pos], line[pos + 1:]
            value = value.replace("'", "").replace('"', '')
            key, value = key.strip(), value.strip()
            #add_keyword(key, value)

        line = f.readline()


    # Column headers

    l = 0
    units = {}
    pipes = []
    names = []

    while True:

        char1 = line[0:2]

        if char1 <> "#|":
            break

        if l==1: # Column names

            line = line.replace('-', ' ').strip()

            # Find all pipe symbols
            for i, c in enumerate(line):
                if c=='|':
                    pipes.append(i)

            # Find all names
            names = line.replace(" ", "").split("|")[1:-1]

        elif l==2: # Data types

            line = line.replace('-', ' ').strip()

            types = dict(zip(names, \
                line.replace(" ", "").split("|")[1:-1]))

        elif l==3: # Units

            units = dict(zip(names, \
                line.replace(" ", "").split("|")[1:-1]))

        else: # Null values

            pass
        
        line = f.readline()
        l = l + 1

    if len(pipes) <> len(names) + 1:
        print "An error occured while reading the c2d table"

    if len(units)==0:
        for name in names:
            units[name]=''

    # if len(nulls)==0:
    #             for name in names:
    #                 nulls[name]=''
    
    # Data

    array = {}
    for name in names:
        array[name] = []


    while True:

        if line.strip() == '':
            break

        for i in range(len(pipes)-1):

            first, last = pipes[i] + 1, pipes[i + 1]

            if definition==1:
                last = last + 1
                if first==1:
                    first=0
            elif definition==2:
                first = first - 1

            if i + 1==len(pipes)-1:
                item = line[first:].strip()
            else:
                item = line[first:last].strip()

            # if item == nulls[names[i]]:
            #                     item = 'NaN'
            array[names[i]].append(item)

        line = f.readline()

    # Convert to numpy arrays and add to table
    for name in names:
        array[name] = np.array(array[name], \
            dtype=type_dict[types[name]])
        atpytbl.add_column(name, array[name], \
            unit=units[name])
