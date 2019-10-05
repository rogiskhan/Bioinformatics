#this script creates a list of zeros for each human chromosome.
#
if __name__=="__main__": 
    #initialize a list with the number of bases for each chromosome (1 to 22, X, Y)
    list=[248956422, 242193529,	198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415]
    #create txt file for each chromosome ('till 22)
    for j in range(22):
        au=j+1
        au.__str__() #convert int into a string to put it in the name
        ch="chr%s.txt" %au #create the name string
        
        f=open("%s" %ch, "w") #open the ch-th file
        #write zeros
        for i in range(list[j]):
            f.write("0")
        f.close
    #same as before but for x and y
    f=open("chrX.txt", "w")
    for i in range(list[22]):
        f.write("0")
    f.close
    f=open("chrY.txt", "w")
    for i in range(list[23]):
        f.write("0")
    f.close
pass