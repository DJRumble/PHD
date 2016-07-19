""" A module for completing the trapizum rule as a method of calculating an integral """ 

def trap(X,Y):
    #Takes a list of X points and Y(x). Both lists must be the same length - bin sizes are uneven

    i = 0
    B = 0
    for i in range(len(X)-1):
        
        a1 = (X[i+1]-X[i])
        a2 = (Y[i+1]+Y[i])
        A = a1*a2
        #print 'MODULE: a1 = ', a1
        #print 'MODULE: ', Y[i+1], ' & ', Y[i]
        #print 'MODULE: a2 = ', a2
        B = B + A
        i = i + 1

    I = 0.5 * B

    return I 

if __name__ == "__main__":
    trap(X,Y)
        
