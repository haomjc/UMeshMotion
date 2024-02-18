import numpy  
my_matrix = numpy.loadtxt(open("E:\\temp\\UMeshMotion\\RoughSurface_135.csv","rb"),delimiter=",",skiprows=0, dtype=float)
print my_matrix

p = mdb.models['Job-230515_1_wear_1_6_M_4_T'].parts['ROLLER']
n = p.nodes

for i in range(31611,44436,135):
    for j in range(i,i+135):

        nodes = n[j-1:j]
        p.editNode(nodes=nodes, offset1=0, offset2=0, offset3=my_matrix[(i-31611)/135][j-i])