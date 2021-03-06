'''
You just create a table object, specifying the number of columns (6), their justifications (justs), the table caption, and \ref label to use.  
Next, add the header labels.  
Then add data, which can be a list of strings (col1), a list of numbers (col2) or a list of number and error pairs (col3).   
The argument sigfigs lets you set the number of significant figures for the data in the table.  
It will round up to this for single numbers and round up the errors and format the numbers appropriately.
'''

import Table
fout = open('mytable.tex','w')
t = Table.Table(6, justs='lrc', caption='Awesome results', label="tab:label")
t.add_header_row(['obj', 'X', '$\\beta$'])
col1 = ['obj1','obj2','obj3']
col2 = [0.001,0.556,10.56]   # just numbers
col3 = [[0.12345,0.1],[0.12345,0.01],[0.12345,0.001]]
t.add_data([col1,col2,col3], sigfigs=2)
t.print_table(fout)
fout.close()
