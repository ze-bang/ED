import matplotlib.pyplot as plt
import numpy as np
from opt_einsum import contract
import sys

z = np.array([[1,1,1],[1,-1,-1],[-1,1,-1], [-1,-1,1]])/np.sqrt(3)
y = np.array([[0,-1,1],[0,1,-1],[0,-1,-1], [0,1,1]])/np.sqrt(2)
x = np.array([[-2,1,1],[-2,-1,-1],[2,1,-1], [2,-1,1]])/np.sqrt(6)

coord = np.array([x,y,z])

pyrochore_unit_cell = np.array([[-1,-1,-1],[-1,1,1],[1,-1,1],[1,1,-1]])/8
basis = np.array([[0,1,1],[1,0,1],[1,1,0]])/2
BasisBZA = np.array([2*np.pi*np.array([-1,1,1]),2*np.pi*np.array([1,-1,1]),2*np.pi*np.array([1,1,-1])])

dim1 = 1
dim2 = 1
dim3 = 1

coord = np.zeros((dim1*dim2*dim3*4,3))

for i in range(dim1):
    for j in range(dim2):
        for k in range(dim3):
            for u in range(4):
                coord[i*dim2*dim3*4+j*dim3*4+k*4+u] = basis[0]*i + basis[1]*j + basis[2]*k + pyrochore_unit_cell[u]


def readTwoBody(filename):
    f = open(filename, 'r')

dir = sys.argv[1]

twobody = np.loadtxt(dir+"/output/zvo_cisajscktalt_eigen0.dat", dtype=np.float64)
onebody = np.loadtxt(dir+"/output/zvo_cisajs_eigen0.dat", dtype=np.float64)
# ind = np.array([1,3,5,7], dtype=int)
# spin_info = f[:,ind]

def gSF(q, v):
    M = np.zeros((len(q),4,4))
    a = np.cross(q, v)
    a = contract('ia, i->ia', a, 1/np.linalg.norm(a, axis=1))
    for i in range(4):
        for j in range(4):
            M[:,i,j] = contract('a, ka, b, kb->k',z[i], a, z[j], a)
    return M

def gNSF(q, v):
    M = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            M[i,j] = contract('a,a, b, b->',z[i], v, z[j], v)
    M = contract('k ,ab->kab', np.ones(len(q)), M)
    return M

def proj(k):
    delta = contract('ar, br->ab', z, z)
    return delta - contract('ar, ir, bp, ip, i->iab', z, k,
                    z, k, 1/contract('ik,ik->i', k ,k))
def projx(k):
    delta = contract('ar, br->ab', x, x)
    return delta - contract('ar, ir, bp, ip, i->iab', x, k,
                    x, k, 1/contract('ik,ik->i', k ,k))
def projy(k):
    delta = contract('ar, br->ab', y, y)
    return delta - contract('ar, ir, bp, ip, i->iab', y, k,
                    y, k, 1/contract('ik,ik->i', k ,k))

def Sz(k, f):
    indpos = np.arange(0, len(f), 4)
    indneg = np.arange(3, len(f), 4)
    CupCup = f[indpos,4] + 1j * f[indpos,5]
    CupCdown = f[indneg,4] + 1j * f[indneg,5]

    Szz = (CupCup - CupCdown)/2

    ind_pos_1 = np.array(f[indpos,0], dtype=int)

    here_coord_1 = coord[ind_pos_1]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1))

    return np.real(contract('n,in->i', Szz, ffact))

def SzSz(k,f):
    Sk = Sz(k, f)
    return contract('i,i->i',Sk, np.conj(Sk))


def Szz(k, f):
    indpos = np.concatenate((np.arange(0, len(f), 8), np.arange(3, len(f), 8)))
    indneg = np.concatenate((np.arange(1, len(f), 8), np.arange(2, len(f), 8)))
    CupCup = f[indpos,8] + 1j * f[indpos,9]
    CupCdown = f[indneg,8] + 1j * f[indneg,9]

    Szz = (CupCup - CupCdown)/4

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return np.real(contract('n,in->i', Szz, ffact))
def Spm(k, f):
    indpos = np.arange(5, len(f), 8)
    Spm_Smp = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in->i', Spm_Smp, ffact)
def Smp(k, f):
    indpos = np.arange(4, len(f), 8)
    Spm_Smp = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in->i', Spm_Smp, ffact)
def Spp(k, f):
    indpos = np.arange(7, len(f), 8)
    Spp_Smm = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in->i', Spp_Smm, ffact)
def Smm(k, f):
    indpos = np.arange(6, len(f), 8)
    Spp_Smm = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in->i', Spp_Smm, ffact)

def Szz_global(k, f):
    indpos = np.concatenate((np.arange(0, len(f), 8), np.arange(3, len(f), 8)))
    indneg = np.concatenate((np.arange(1, len(f), 8), np.arange(2, len(f), 8)))
    CupCup = f[indpos,8] + 1j * f[indpos,9]
    CupCdown = f[indneg,8] + 1j * f[indneg,9]

    Szz = (CupCup - CupCdown)/4

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    P = proj(k)
    sublat_ind_1 = np.mod(ind_pos_1, 4)
    sublat_ind_2 = np.mod(ind_pos_2, 4)

    realProj = np.zeros((len(k),len(sublat_ind_1)))

    for i in range(len(sublat_ind_1)):
        realProj[:, i] = P[:,sublat_ind_1[i], sublat_ind_2[i]]

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return np.real(contract('n,in, in->i', Szz, ffact, realProj))

def Szz_NSF(k, f, v):
    indpos = np.concatenate((np.arange(0, len(f), 8), np.arange(3, len(f), 8)))
    indneg = np.concatenate((np.arange(1, len(f), 8), np.arange(2, len(f), 8)))
    CupCup = f[indpos,8] + 1j * f[indpos,9]
    CupCdown = f[indneg,8] + 1j * f[indneg,9]

    Szz = (CupCup - CupCdown)/4

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    P = gNSF(k, v)
    sublat_ind_1 = np.mod(ind_pos_1, 4)
    sublat_ind_2 = np.mod(ind_pos_2, 4)

    realProj = np.zeros((len(k),len(sublat_ind_1)))

    for i in range(len(sublat_ind_1)):
        realProj[:,i] = P[:,sublat_ind_1[i], sublat_ind_2[i]]

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return np.real(contract('n,in, in->i', Szz, ffact, realProj))


def Spm_global(k, f,proj_method):
    indpos = np.arange(5, len(f), 8)
    Spm_Smp = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    P = proj_method(k)
    sublat_ind_1 = np.mod(ind_pos_1, 4)
    sublat_ind_2 = np.mod(ind_pos_2, 4)

    realProj = np.zeros((len(k),len(sublat_ind_1)))

    for i in range(len(sublat_ind_1)):
        realProj[:, i] = P[:,sublat_ind_1[i], sublat_ind_2[i]]

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in,in->i', Spm_Smp, ffact, realProj)
def Smp_global(k, f, proj_method):
    indpos = np.arange(4, len(f), 8)
    Spm_Smp = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)
    P = proj_method(k)
    sublat_ind_1 = np.mod(ind_pos_1, 4)
    sublat_ind_2 = np.mod(ind_pos_2, 4)

    realProj = np.zeros((len(k),len(sublat_ind_1)))

    for i in range(len(sublat_ind_1)):
        realProj[:, i] = P[:,sublat_ind_1[i], sublat_ind_2[i]]
    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in,in->i', Spm_Smp, ffact, realProj)
def Spp_global(k, f, proj_method):
    indpos = np.arange(7, len(f), 8)
    Spp_Smm = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    P = proj_method(k)
    sublat_ind_1 = np.mod(ind_pos_1, 4)
    sublat_ind_2 = np.mod(ind_pos_2, 4)

    realProj = np.zeros((len(k),len(sublat_ind_1)))

    for i in range(len(sublat_ind_1)):
        realProj[:, i] = P[:,sublat_ind_1[i], sublat_ind_2[i]]

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in, in->i', Spp_Smm, ffact, realProj)
def Smm_global(k, f, proj_method):
    indpos = np.arange(6, len(f), 8)
    Spp_Smm = f[indpos,8] + 1j * f[indpos,9]

    ind_pos_1 = np.array(f[indpos,0], dtype=int)
    ind_pos_2 = np.array(f[indpos,4], dtype=int)

    P = proj_method(k)
    sublat_ind_1 = np.mod(ind_pos_1, 4)
    sublat_ind_2 = np.mod(ind_pos_2, 4)

    realProj = np.zeros((len(k),len(sublat_ind_1)))

    for i in range(len(sublat_ind_1)):
        realProj[:, i] = P[:,sublat_ind_1[i], sublat_ind_2[i]]

    here_coord_1 = coord[ind_pos_1]
    here_coord_2 = coord[ind_pos_2]

    ffact = np.exp(1j*contract('ik, nk->in', k, here_coord_1-here_coord_2))

    return contract('n,in, in->i', Spp_Smm, ffact, realProj)
def Sxx(k, f):
    return np.real(Spm(k,f) + Smp(k,f) + Spp(k, f) + Smm(k, f))/4
def Syy(k, f):
    return np.real(Spm(k,f) + Smp(k,f) - Spp(k, f) - Smm(k, f))/4
def Sxy(k, f):
    return np.real((-Spm(k,f) + Smp(k,f) + Spp(k, f) - Smm(k, f))/(4j))
def Syx(k, f):
    return np.real((Spm(k,f) - Smp(k,f) + Spp(k, f) - Smm(k, f))/(4j))
def Sxx_global(k, f):
    return np.real(Spm_global(k,f, projx) + Smp_global(k,f, projx) + Spp_global(k, f, projx) + Smm_global(k, f, projx))/4
def Syy_global(k, f):
    return np.real(Spm_global(k,f, projy) + Smp_global(k,f, projy) - Spp_global(k, f, projy) - Smm_global(k, f, projy))/4
def hhltoK(H, L):
    return np.einsum('ij,k->ijk',H, 2*np.array([np.pi,np.pi,0])) \
        + np.einsum('ij,k->ijk',L, 2*np.array([0.,0.,np.pi]))

def hk2ktoK(H, L):
    return np.einsum('ij,k->ijk',H, 2*np.array([np.pi,-np.pi,0])) \
        + np.einsum('ij,k->ijk',L, 2*np.array([np.pi,np.pi,-2*np.pi]))


def genHHLK(Hr, Lr, nK):
    H = np.linspace(-Hr, Hr, nK)
    L = np.linspace(-Lr, Lr, nK)
    A, B = np.meshgrid(H, L)
    return hhltoK(A, B).reshape((nK * nK, 3))
def genHnH2KK(Hr, Lr, nK):
    H = np.linspace(-Hr, Hr, nK)
    L = np.linspace(-Lr, Lr, nK)
    A, B = np.meshgrid(H, L)
    return hk2ktoK(A, B).reshape((nK * nK, 3))


nK = 100

if sys.argv[2] == "111":
    Hr = 3
    Lr = 3
    K = genHnH2KK(Hr, Lr, nK)
else:
    Hr = 2.5
    Lr = 2.5
    K = genHHLK(Hr, Lr, nK)


SSSF = Szz(K, twobody)

SSSF = SSSF.reshape((nK,nK))

C = plt.imshow(SSSF, interpolation="lanczos", origin='lower', extent=[-Hr, Hr, -Lr, Lr], aspect='equal')
# plt.colorbar(C)
# plt.ylabel(r'$(0,0,L)$')
# plt.xlabel(r'$(H,H,0)$')
plt.xlim([-Hr, Hr])
plt.ylim([-Lr, Lr])
plt.savefig(dir+"/Szz.png")
plt.clf()

SSSF = Szz_global(K, twobody)

SSSF = SSSF.reshape((nK,nK))

C = plt.imshow(SSSF, interpolation="lanczos", origin='lower', extent=[-Hr, Hr, -Lr, Lr], aspect='equal')
# plt.colorbar(C)
# plt.ylabel(r'$(K,K,-2K)$')
# plt.xlabel(r'$(H,-H,0)$')
plt.xlim([-Hr, Hr])
plt.ylim([-Lr, Lr])
plt.savefig(dir+"/Szz_global.png")
plt.clf()


SSSF = Szz_NSF(K, twobody, np.array([1,-1,0])/np.sqrt(2))

SSSF = SSSF.reshape((nK,nK))

C = plt.imshow(SSSF, interpolation="lanczos", origin='lower', extent=[-Hr, Hr, -Lr, Lr], aspect='equal')
plt.colorbar(C)
# plt.ylabel(r'$(K,K,-2K)$')
# plt.xlabel(r'$(H,-H,0)$')
plt.xlim([-Hr, Hr])
plt.ylim([-Lr, Lr])
plt.savefig(dir+"/Szz_NSF.png")
plt.clf()