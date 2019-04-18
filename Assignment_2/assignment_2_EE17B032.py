import sys

import numpy as np
from numpy import cos,sin

class one_port_element():
    ''' General class for all the one port elements '''
    def __init__(self,line):
        self.line = line
        self.tokens = self.line.split()
        self.name = element_type(self.tokens[0])
        self.from_node = self.tokens[1]
        self.to_node = self.tokens[2]
        if len(self.tokens) == 5:
            self.type = 'dc'
            self.value = float(self.tokens[4])
        elif len(self.tokens) == 6:
            self.type = 'ac'
            Vm = float(self.tokens[4])/2
            phase = float(self.tokens[5])
            real = Vm*cos(phase)
            imag = Vm*sin(phase)
            self.value = complex(real,imag)
        else:
            self.type = 'RLC'
            self.value = float(self.tokens[3])

def remove_comments(cir_def):
    ''' Removes the comments in the circuit definition lines '''
    mod_def = []
    for line in cir_def:
        ind = len(line)-1
        for i in range(len(line)):
            if line[i] == '#' or line[i] =='#\n' :
                ind = i
                break
        mod_def.append(line[:ind])
    return mod_def

def element_type(token):
    ''' Gets the element name '''
    if token[0] == 'R':
        return 'resistor'
    elif token[0] == 'L':
        return 'inductor'
    elif token[0] == 'C':
        return 'capacitor'
    elif token[0] == 'V':
        return 'ind voltage source'
    elif token[0] == 'I':
        return 'ind current source'

def convert_to_dict(cir_def):
    ''' Returns a dictionary of nodes from the circuit definition. By default, the GND node is assigned value 0. '''
    diction = {}
    nodes = [one_port_element(line).from_node for line in cir_def]
    nodes.extend([one_port_element(line).to_node for line in cir_def])
    ind = 1
    nodes = list(set(nodes))
    for node in nodes:
        if node == 'GND' :
            diction[node] = 0
        else :
            diction[node] = ind
            ind += 1
    return diction

def get_key(diction,value):
    ''' Gets the corresponding key for a value in the dictionary '''
    for key in diction.keys():
        if diction[key] == value :
            return key
def make_dict(cir_def,element):
    ''' Makes a dictionary for each component of the particular type of element ''' 
    e = element
    volt_dict = {}
    volt_names = [one_port_element(line).tokens[0] for line in cir_def if one_port_element(line).tokens[0][0].lower()== e]
    for ind,name in enumerate(volt_names):
        volt_dict[name] = ind
    return volt_dict


AC = '.ac' 
def get_freq(lines):
    ''' Returns the frequency of the source '''
    w = 0
    for line in lines:
        if line[:len(AC)] == '.ac' :
            w = float(line.split()[2])
    return w


def matrix_dims(cir_def):
    ''' Returns a tuple : number of nodes, number of voltage sources '''
    volt_ind = [i for i in range(len(cir_def)) if cir_def[i].split()[0][0] == 'V']
    k = len(volt_ind)
    n = len(convert_to_dict(cir_def))
    return n,k
    
def node_finder(cir_def,node_key,diction):
    ''' Finds the lines and position ie from/to of the given node '''
    inds = [(i,j) for i in range(len(cir_def)) for j in range(len(cir_def[i].split())) if cir_def[i].split()[j] in diction.keys() if diction[cir_def[i].split()[j]] == node_key]
    return inds


def mod_matrix(cir_def,w,node_key,diction,volt_dict,ind_dict,M,b):
    ''' Updates the M and b matrix for the given node '''
    inds = node_finder(cir_def,node_key,diction)
    n,k = matrix_dims(cir_def)
    for ind in inds:
        #getting all the attributes of the element using the class definition
        element = one_port_element(cir_def[ind[0]])
        element_name = cir_def[ind[0]].split()[0]
        if element_name[0] == 'R':
            if ind[1] == 1:
                #node is the from_node
                adj_key = diction[element.to_node]
                M[node_key,node_key] += 1/(element.value)
                M[node_key,adj_key] -= 1/(element.value)
                    
            if ind[1] == 2 :
                # node is the to_node
                adj_key = diction[element.from_node]
                M[node_key,node_key] += 1/(element.value)
                M[node_key,adj_key] -= 1/(element.value)
                
        if element_name[0] == 'C' :
            if ind[1]== 1:
                # node is the from_node
                adj_key = diction[element.to_node]
                M[node_key,node_key] += complex(0, 2*np.pi*w*(element.value))
                M[node_key,adj_key] -= complex(0, 2*np.pi*w*(element.value))
            if ind[1] == 2 :
                # node is the to_node
                adj_key = diction[element.from_node]
                M[node_key,node_key] += complex(0, 2*np.pi*w*(element.value))
                M[node_key,adj_key] -= complex(0, 2*np.pi*w*(element.value))

        if element_name[0] == 'L' :
            try:
                if ind[1]== 1:
                    adj_key = diction[element.to_node]
                    M[node_key,node_key] -= complex(0,1/(2*np.pi*w*element.value))
                    M[node_key,adj_key] += complex(0,1/(2*np.pi*w*element.value))
                if ind[1] == 2 :
                    adj_key = diction[element.from_node]
                    M[node_key,node_key] -= complex(0,1/(2*np.pi*w*element.value))
                    M[node_key,adj_key] += complex(0,1/(2*np.pi*w*element.value))
            except ZeroDivisionError:
                index = ind_dict[element_name]
                if ind[1]== 1:
                    adj_key = diction[element.to_node]
                    M[node_key,n+k+index] += 1 
                    M[n+k+index,node_key] -= 1
                    b[n+k+index] = 0
                if ind[1]== 2:
                    M[node_key,n+k+index] -= 1
                    M[n+k+index,node_key] += 1
                    b[n+k+index] = 0
        if element_name[0] == 'V' :
            index = volt_dict[element_name]
            if ind[1]== 1:
                adj_key = diction[element.to_node]
                M[node_key,n+index] += 1
                M[n+index,node_key] -= 1
                b[n+index] = element.value
            if ind[1] == 2 :
                adj_key = diction[element.from_node]
                M[node_key,n+index] -= 1
                M[n+index,node_key] +=1
                b[n+index] = element.value

        if element_name[0] == 'I' :
            if ind[1]== 1:
                b[node_key] -= element.value
            if ind[1] == 2 :
                b[node_key] += element.value

''' Main function '''
try:
    if( len(sys.argv) != 2):
        sys.exit('Invalid number of arguments')
    name = sys.argv[1]
    f = open(name)
    lines = f.readlines()
    f.close()
    w =get_freq(lines)
    # getting the circuit definition block :
    for ind,string in enumerate(lines):
        try : 
            if string.split()[0] == '.circuit\n' or string.split()[0] == '.circuit':
                start_ind = ind
            elif string.split()[0] == '.end' or string.split()[0]  == '.end\n':
                end_ind = ind
        except IndexError: 
            continue
    cir_def = lines[start_ind+1:end_ind]
    cir_def = remove_comments(cir_def)
    diction = convert_to_dict(cir_def)
    volt_dict = make_dict(cir_def,'v')
    ind_dict = make_dict(cir_def,'l')
    n,k = matrix_dims(cir_def)
    dim = n+k
    M = np.zeros((dim,dim),dtype=np.complex)
    b = np.zeros(dim,dtype=np.complex)
    dc_flag = False 
    if w == 0:
        dc_flag = True
        M = np.zeros((dim+len(ind_dict),dim+len(ind_dict)),dtype=np.complex)
        b = np.zeros(dim+len(ind_dict),dtype=np.complex)

    for i in range(len(diction)):
        mod_matrix(cir_def,w,i,diction,volt_dict,ind_dict,M,b)
    M[0] = 0
    M[0,0] =1
    print('The M matrix is :\n',M)
    print('The b matrix is :\n',b)
    print('The node dictionary is :',diction)
    try:
        x = np.linalg.solve(M,b)    
    except Exception:
        print('The incidence matrix cannot be inverted as it is singular. Please provide a valid circuit definition')
        sys.exit()
    for i in range(n):
        print("The voltage at node {} is {}".format(get_key(diction,i),x[i]))
    for j in range(k):
        print('The current through source {} is {}'.format(get_key(volt_dict,j),x[n+j]))
    if dc_flag:
        for i in range(len(ind_dict)):
            print("The current through inductor {} is {}".format(get_key(ind_dict,i),x[n+k+i]))
    print('Voltage convention : From node of the voltage source is at a lower potential')     

except Exception:
   print('Invalid file. Please make sure that the circuit definition block is well defined and all component value are in scientific notation.')
