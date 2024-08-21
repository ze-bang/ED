import numpy as np

z = np.array([np.array([1,1,1])/np.sqrt(3), np.array([1,-1,-1])/np.sqrt(3), np.array([-1,1,-1])/np.sqrt(3), np.array([-1,-1,1])/np.sqrt(3)])

def indices(i,j,k,u,d1, d2, d3):
    if u == 0:
        return np.array([[i, j, k, 1], [i, j, k, 2], [i, j, k, 3], [np.mod(i-1, d1), j, k, 1], [i, np.mod(j-1, d2), k, 2], [i, j, np.mod(k-1, d3), 3]])
    elif u == 1:
        return np.array([[i, j, k, 0], [i, j, k, 2], [i, j, k, 3], [np.mod(i+1, d1), j, k, 0],[np.mod(i+1, d1), np.mod(j-1, d2), k, 2], [np.mod(i+1, d1), j, np.mod(k-1, d3), 3]])
    elif u == 2:
        return np.array([[i, j, k, 0], [i, j, k, 1], [i, j, k, 3], [i, np.mod(j+1, d2), k, 0], [np.mod(i-1, d1), np.mod(j+1, d2), k, 1], [i, np.mod(j+1, d2), np.mod(k-1, d3), 3]])
    elif u == 3:
        return np.array([[i, j, k, 0], [i, j, k, 1], [i, j, k, 2], [i, j, np.mod(k+1, d3), 0], [np.mod(i-1, d1), j, np.mod(k+1, d3), 1], [i, np.mod(j-1, d2), np.mod(k+1, d3), 2]])


dim1 = 1
dim2 = 2
dim3 = 2


con = np.zeros((dim1, dim2, dim3, 4))

Jxx, Jyy, Jzz = 0.6, 1.0, 0.6
Jpm = -(Jxx+Jyy)/4
Jpmpm = (Jxx-Jyy)/4
h = 0.0*2.24
fielddir = np.array([1,1,0])/np.sqrt(2)
B = np.einsum('r, ir->i', h*fielddir,z)


def flattenIndex(Indx):
    temp = np.zeros(6)
    for i in range(6):
        temp[i] = Indx[i][0]*dim2*dim3*4 + Indx[i][1]*dim3*4 + Indx[i][2]*4 + Indx[i][3]
    return temp

def genNN_list(d1,d2,d3):
    NN_list = np.zeros((d1*d2*d3*4, 6))
    for i in range(d1):
        for j in range(d2):
            for k in range(d3):
                for u in range(4):
                    NN_list[i*d2*d3*4+j*d3*4+k*4+u] = flattenIndex(indices(i,j,k,u,d1,d2,d3))
    return NN_list

look_up_table = genNN_list(dim1, dim2, dim3)

def HeisenbergNN(Jzz, Jpm, Jpmpm, indx1, indx2):
    return np.array([[indx1, 0, indx1, 0, indx2, 0, indx2, 0, Jzz/8, 0],
                     [indx1, 0, indx1, 0, indx2, 1, indx2, 1, -Jzz/8, 0],
                     [indx1, 1, indx1, 1, indx2, 0, indx2, 0, -Jzz/8, 0],
                     [indx1, 1, indx1, 1, indx2, 1, indx2, 1, Jzz/8, 0],

                     [indx1, 1, indx1, 0, indx2, 0, indx2, 1, -Jpm/2, 0],
                     [indx1, 0, indx1, 1, indx2, 1, indx2, 0, -Jpm/2, 0],

                     [indx1, 1, indx1, 0, indx2, 1, indx2, 0, Jpmpm/2, 0],
                     [indx1, 0, indx1, 1, indx2, 0, indx2, 1, Jpmpm/2, 0]])

def Zeeman(h, indx):
    here = h[indx % 4]
    return np.array([[indx, 0, indx, 0, here/2, 0],
                     [indx, 1, indx, 1, -here/2, 0]])   


interALL = []
transfer = []


for i in range(len(look_up_table)):
    transfer.append(Zeeman(B, i))
    for j in range(len(look_up_table[i])):
        interALL.append(HeisenbergNN(Jzz, Jpm, Jpmpm, i, look_up_table[i][j]))

interALL = np.array(interALL).reshape(-1, 10)
transfer = np.array(transfer).reshape(-1, 6)

def write_interALL(interALL, file_name):
    num_param = len(interALL)
    f        = open(file_name, 'wt')
    f.write("==================="+"\n")
    f.write("num "+"{0:8d}".format(num_param)+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    for i in range(num_param):
        f.write(" {0:8d} ".format(int(interALL[i,0])) \
        +" {0:8d}   ".format(int(interALL[i,1]))     \
        +" {0:8d}   ".format(int(interALL[i,2]))     \
        +" {0:8d}   ".format(int(interALL[i,3]))     \
        +" {0:8d}   ".format(int(interALL[i,4]))     \
        +" {0:8d}   ".format(int(interALL[i,5]))     \
        +" {0:8d}   ".format(int(interALL[i,6]))     \
        +" {0:8d}   ".format(int(interALL[i,7]))     \
        +" {0:8f}   ".format(interALL[i,8]) \
        +" {0:8f}   ".format(interALL[i,9]) \
        +"\n")
    f.close()

def write_transfer(interALL, file_name):
    num_param = len(interALL)
    f        = open(file_name, 'wt')
    f.write("==================="+"\n")
    f.write("num "+"{0:8d}".format(num_param)+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    f.write("==================="+"\n")
    for i in range(num_param):
        f.write(" {0:8d} ".format(int(interALL[i,0])) \
        +" {0:8d}   ".format(int(interALL[i,1]))     \
        +" {0:8d}   ".format(int(interALL[i,2]))     \
        +" {0:8d}   ".format(int(interALL[i,3]))     \
        +" {0:8f}   ".format(interALL[i,4])     \
        +" {0:8f}   ".format(interALL[i,5])
        +"\n")
    f.close()


output_dir = "./"
max_site = dim1*dim2*dim3*4
All_N = max_site
exct = 1


write_interALL(interALL, 'InterAll.def')
write_transfer(transfer, 'Trans.def')

f        = open("{}".format(output_dir)+"calcmod_cg.def", 'wt')
f.write("  #CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution"+"\n")
f.write("  #CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC"+"\n")
f.write("  #Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart"+"\n")
f.write("  #CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save"+"\n")
f.write("  CalcType   3"+"\n")
f.write("  CalcModel   4"+"\n")
f.write("  ReStart   0"+"\n")
f.write("  CalcSpec   0"+"\n")
f.write("  CalcEigenVec   0"+"\n")
f.write("  InitialVecType   0"+"\n")
f.write("  InputEigenVec   0"+"\n")
f.write("  OutputEigenVec   0"+"\n")
f.close()

f        = open("{}".format(output_dir)+"calcmod_tpq.def", 'wt')
f.write("  #CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution"+"\n")
f.write("  #CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC"+"\n")
f.write("  #Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart"+"\n")
f.write("  #CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save"+"\n")
f.write("  CalcType   1"+"\n")
f.write("  CalcModel   4"+"\n")
f.write("  ReStart   0"+"\n")
f.write("  CalcSpec   0"+"\n")
f.write("  CalcEigenVec   0"+"\n")
f.write("  InitialVecType   0"+"\n")
f.write("  InputEigenVec   0"+"\n")
f.write("  OutputEigenVec   0"+"\n")
f.close()


f        = open("{}".format(output_dir)+"namelist_cg.def", 'wt')
f.write("  ModPara       modpara.def"+"\n")
f.write("  CalcMod       calcmod_cg.def"+"\n")
f.write("  LocSpin       locspn.def"+"\n")
f.write("  Trans         Trans.def"+"\n")
f.write("  InterAll      InterAll.def"+"\n")
f.write("  OneBodyG      greenone.def"+"\n")
f.write("  TwoBodyG      greentwo.def"+"\n")
f.close()

f        = open("{}".format(output_dir)+"namelist_tpq.def", 'wt')
f.write("  ModPara       modpara.def"+"\n")
f.write("  CalcMod       calcmod_tpq.def"+"\n")
f.write("  LocSpin       locspn.def"+"\n")
f.write("  Trans         Trans.def"+"\n")
f.write("  InterAll      InterAll.def"+"\n")
f.write("  OneBodyG      greenone.def"+"\n")
f.write("  TwoBodyG      greentwo.def"+"\n")
f.close()
 

f        = open("./{}/".format(output_dir)+"modpara.def", 'wt')
f.write("--------------------  "+"\n")
f.write("Model_Parameters   0  "+"\n")
f.write("--------------------  "+"\n")
f.write("HPhi_Cal_Parameters  "+"\n")
f.write("--------------------  "+"\n")
f.write("CDataFileHead  zvo  "+"\n")
f.write("CParaFileHead  zqp  "+"\n")
f.write("--------------------  "+"\n")
f.write("Nsite          {}".format(max_site)+"\n")
f.write("Ncond          {}".format(max_site)+"\n")
f.write("Lanczos_max    2000   "+"\n")
f.write("initial_iv     -1    "+"\n")
f.write("exct           {}".format(exct)+"\n")
f.write("LanczosEps     14   "+"\n")
f.write("LanczosTarget  2   "+"\n")
f.write("LargeValue     30  "+"\n")
f.write("NumAve         5   "+"\n")
f.write("ExpecInterval  20  "+"\n")
f.close()

num_loc  = max_site
f        = open("{}".format(output_dir)+"locspn.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_loc)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
for all_i in range(0,All_N):
     f.write(" {0:8d} ".format(all_i)+" 1 " \
     +"\n")
f.close()

num_green  = 4*All_N
f        = open("{}".format(output_dir)+"greenone.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_green)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
for all_i in range(0,All_N):
    f.write(" {0:8d} ".format(all_i)+" 0 " \
        +" {0:8d} ".format(all_i)+" 0 "     \
        +"\n")
    f.write(" {0:8d} ".format(all_i)+" 0 " \
        +" {0:8d} ".format(all_i)+" 1 "     \
        +"\n")
    f.write(" {0:8d} ".format(all_i)+" 1 " \
        +" {0:8d} ".format(all_i)+" 0 "     \
        +"\n")
    f.write(" {0:8d} ".format(all_i)+" 1 " \
        +" {0:8d} ".format(all_i)+" 1 "     \
        +"\n")
f.close()

num_green  = 8*All_N*All_N
f        = open("{}".format(output_dir)+"greentwo.def", 'wt')
f.write("==================="+"\n")
f.write("loc "+"{0:8d}".format(num_green)+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
f.write("==================="+"\n")
for all_i in range(0,All_N):
    for all_j in range(0,All_N):
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 0 " \
           +" {0:8d} ".format(all_i)+" 1 "     \
           +" {0:8d} ".format(all_j)+" 0 "     \
           +" {0:8d}   ".format(all_j)+" 1 "   \
           +"\n")
        f.write(" {0:8d} ".format(all_i)+" 1 " \
           +" {0:8d} ".format(all_i)+" 0 "     \
           +" {0:8d} ".format(all_j)+" 1 "     \
           +" {0:8d}   ".format(all_j)+" 0 "   \
           +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 1 " \
        #    +" {0:8d} ".format(all_i)+" 0 "     \
        #    +" {0:8d} ".format(all_j)+" 0 "     \
        #    +" {0:8d}   ".format(all_j)+" 0 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 0 " \
        #    +" {0:8d} ".format(all_i)+" 1 "     \
        #    +" {0:8d} ".format(all_j)+" 0 "     \
        #    +" {0:8d}   ".format(all_j)+" 0 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 0 " \
        #    +" {0:8d} ".format(all_i)+" 0 "     \
        #    +" {0:8d} ".format(all_j)+" 1 "     \
        #    +" {0:8d}   ".format(all_j)+" 0 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 0 " \
        #    +" {0:8d} ".format(all_i)+" 0 "     \
        #    +" {0:8d} ".format(all_j)+" 0 "     \
        #    +" {0:8d}   ".format(all_j)+" 1 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 0 " \
        #    +" {0:8d} ".format(all_i)+" 1 "     \
        #    +" {0:8d} ".format(all_j)+" 1 "     \
        #    +" {0:8d}   ".format(all_j)+" 1 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 1 " \
        #    +" {0:8d} ".format(all_i)+" 0 "     \
        #    +" {0:8d} ".format(all_j)+" 1 "     \
        #    +" {0:8d}   ".format(all_j)+" 1 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 1 " \
        #    +" {0:8d} ".format(all_i)+" 1 "     \
        #    +" {0:8d} ".format(all_j)+" 0 "     \
        #    +" {0:8d}   ".format(all_j)+" 1 "   \
        #    +"\n")
        # f.write(" {0:8d} ".format(all_i)+" 1 " \
        #    +" {0:8d} ".format(all_i)+" 1 "     \
        #    +" {0:8d} ".format(all_j)+" 1 "     \
        #    +" {0:8d}   ".format(all_j)+" 0 "   \
        #    +"\n")
f.close()