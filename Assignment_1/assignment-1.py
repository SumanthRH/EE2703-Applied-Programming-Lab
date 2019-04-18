import sys
try :
    if( len(sys.argv) != 2):
        sys.exit('Invalid number of arguments')
    name = sys.argv[1]
    f = open(name)
    lines = f.readlines()
    f.close()
    # getting the cricuit definition block :
    for ind,string in enumerate(lines):
        try : 
            if string.split()[0] == '.circuit\n' or string.split()[0] == '.circuit':
                start_ind = ind
            elif string.split()[0] == '.end' or string.split()[0]  == '.end\n':
                end_ind = ind
        except IndexError: 
            continue
    cir_def = lines[start_ind+1:end_ind]
    # printing in reverse order along with the tokens in reverse :
    print('The circuit definition in reverse order is :')   
    for line in reversed(cir_def) :
        # extracting the tokens :
        tokens = line.split() 
        ind = len(tokens)
        reverse_string = ''
        # ensuring comments are neglected :
        for i in range(len(tokens)):
            if tokens[i] == '#' or tokens[i] =='#\n' :
                ind = i
                break
        for i in range(ind):
            reverse = tokens[ind-(i+1)]
            reverse_string = reverse_string + reverse + ' '        
        print(reverse_string)
except Exception:
    print('Invalid file')


'''if(not(binary & 1)) :
            countz += 1
            maxz = max(maxz,countz)
            maxo = max(maxo,counto) 
            counto = 0
            binary >>= 1

        else:
            counto +=1
            maxo =  max(maxo,counto)
            maxz = max(maxz,countz) 
            countz = 0
            binary >>= 1'''
  