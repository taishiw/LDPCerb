import numpy as np
import math
from scipy import sparse 
from scipy.sparse import csr_matrix
#import cv2
from time import time
#from IPython.display import display, HTML
import scipy
f=open('inputtest.txt','w')
g=open('outputsoftalpha.txt','w')
h=open('outputsoftbeta.txt','w')
o=open('v_sent.txt','w')
p=open('y.txt','w')
numberrow=open('arrayrow.txt','w')
numbercolumn=open('arraycolumn.txt','w')
numberlambda=open('arraylambda.txt','w')
def RegularH(n,d_v,d_c):

    if  n%d_c:
        raise ValueError('d_c must divide n. Help(RegularH) for more info.')

    if d_c <= d_v: 
        raise ValueError('d_c must be greater than d_v. Help(RegularH) for more info.')

    m = (n*d_v)// d_c
    a=m//d_v
    Set=np.zeros((a,n),dtype=int)  
    

    # Filling the first set with consecutive ones in each row of the set 

    for i in range(a):     
        for j in range(i*d_c,(i+1)*d_c): 
            Set[i,j]=1

    #Create list of Sets and append the first reference set
    Sets=[]
    Sets.append(Set.tolist())

    #Create remaining sets by permutations of the first set's columns: 
    i=1
    for i in range(1,d_v):
        newSet = np.transpose(np.random.permutation(np.transpose(Set))).tolist()
        Sets.append(newSet)

    #Returns concatenated list of sest:
    H = np.concatenate(Sets)
    return H

def BinaryProduct(X,Y):
        
    """ Binary Matrices or Matrix-vector product in Z/2Z. Works with scipy.sparse.csr_matrix matrices X,Y too."""
 
    A = X.dot(Y)
    
    if type(A)!=scipy.sparse.csr_matrix:
        return A%2
    
    return A.toarray()%2
def GaussJordan(MATRIX,change=0):

    """ 
    Description:

    Performs the row reduced echelon form of MATRIX and returns it.

    If change = 1, all changes in the MATRIX's rows are applied to identity matrix P: 

    Let A be our parameter MATRIX. refA the reduced echelon form of A. P is the square invertible matrix:

    P.A = Aref.

    -------------------------------------------------
    Parameters: 

    MATRIX: 2D-Array. 
    change : boolean (default = 0)

    ------------------------------------------------

    change = 0  (default)
     >>> Returns 2D-Array Row Reduced Echelon form of Matrix

    change = 1 
    >>> Returns Tuple of 2D-arrays (refMATRIX, P) where P is described above.

    """

    A = np.copy(MATRIX)
    m,n = A.shape
    

    if change:
        P = np.identity(m).astype(int)

    pivot_old = -1 
    for j in range(n):            

        filtre_down = A[pivot_old+1:m,j]
        pivot = np.argmax(filtre_down)+pivot_old+1
        

        if A[pivot,j]:
            pivot_old+=1 
            if pivot_old != pivot:
                aux = np.copy(A[pivot,:])
                A[pivot,:] = A[pivot_old,:]
                A[pivot_old,:] = aux
                if change:
                    aux = np.copy(P[pivot,:])
                    P[pivot,:] = P[pivot_old,:]
                    P[pivot_old,:] = aux

            

            for i in range(m):
                if i!=pivot_old and A[i,j]:
                    if change:
                        P[i,:] = abs(P[i,:]-P[pivot_old,:])
                    A[i,:] = abs(A[i,:]-A[pivot_old,:])


        if pivot_old == m-1:
            break

 
    if change:    
        return A,P 
    return A 

def CodingMatrix(MATRIX,use_sparse=1):

    """ 
    CAUTION: RETURNS tG TRANSPOSED CODING MATRIX. 
    
    Function Applies GaussJordan Algorithm on Columns and rows of MATRIX in order
    to permute Basis Change matrix using Matrix Equivalence.

    Let A be the treated Matrix. refAref the double row reduced echelon Matrix.

    refAref has the form:

    (e.g) : |1 0 0 0 0 0 ... 0 0 0 0|  
            |0 1 0 0 0 0 ... 0 0 0 0| 
            |0 0 0 0 0 0 ... 0 0 0 0| 
            |0 0 0 1 0 0 ... 0 0 0 0| 
            |0 0 0 0 0 0 ... 0 0 0 0| 
            |0 0 0 0 0 0 ... 0 0 0 0| 


    First, let P1 Q1 invertible matrices: P1.A.Q1 = refAref

    We would like to calculate:
    P,Q are the square invertible matrices of the appropriate size so that:

    P.A.Q = J.  Where J is the matrix of the form (having MATRIX's shape):

    | I_p O | where p is MATRIX's rank and I_p Identity matrix of size p.
    | 0   0 |

    Therfore, we perform permuations of rows and columns in refAref (same changes
    are applied to Q1 in order to get final Q matrix)


    NOTE: P IS NOT RETURNED BECAUSE WE DO NOT NEED IT TO SOLVE H.G' = 0 
    P IS INVERTIBLE, WE GET SIMPLY RID OF IT. 
    
    Then
    
    solves: inv(P).J.inv(Q).G' = 0 (1) where inv(P) = P^(-1) and 
    P.H.Q = J. Help(PJQ) for more info.
    
    Let Y = inv(Q).G', equation becomes J.Y = 0 (2) whilst:
    
    J = | I_p O | where p is H's rank and I_p Identity matrix of size p.
        | 0   0 |
    
    Knowing that G must have full rank, a solution of (2) is Y = |  0  | Where k = n-p. 
                                                                 | I-k |
    
    Because of rank-nullity theorem. 
    
    -----------------
    parameters:
    
    H: Parity check matrix. 
    use_sparse: (optional, default True): use scipy.sparse format to speed up calculations
    ---------------
    returns:
    
    tG: Transposed Coding Matrix. 
    
    """


    H = np.copy(MATRIX)
    m,n = H.shape

    if m > n: 
        raise ValueError('MATRIX must have more rows than columns (a parity check matrix)')
    
    if n > 500 and use_sparse:
        sparse = 1
    
    else:
        sparse = 0
    ##### DOUBLE GAUSS-JORDAN:

    Href_colonnes,tQ = GaussJordan(np.transpose(H),1)

    Href_diag = GaussJordan(np.transpose(Href_colonnes))    

    Q=np.transpose(tQ)
    
    k = n - sum(Href_diag.reshape(m*n))

    
    Y = np.zeros(shape=(n,k)).astype(int)
    Y[n-k:,:] = np.identity(k)
    
    if sparse:
        Q = csr_matrix(Q)
        Y = csr_matrix(Y)

    tG = BinaryProduct(Q,Y)

    return tG

def CodingMatrix_systematic(MATRIX,use_sparse = 1):

    """ 
    Description:

    Solves H.G' = 0 and finds the coding matrix G in the systematic form : [I_k  A] by applying permutations on MATRIX.
    
    CAUTION: RETURNS TUPLE (Hp,tGS) WHERE Hp IS A MODIFIED VERSION OF THE GIVEN PARITY CHECK MATRIX, tGS THE TRANSPOSED 
    SYSTEMATIC CODING MATRIX ASSOCIATED TO Hp. YOU MUST USE THE RETURNED TUPLE IN CODING AND DECODING, RATHER THAN THE UNCHANGED 
    PARITY-CHECK MATRIX H. 

    -------------------------------------------------
    Parameters: 

    MATRIX: 2D-Array. Parity-check matrix.
    use_sparse: (optional, default True): use scipy.sparse matrices to speed up calculations if n>100.

    ------------------------------------------------

    >>> Returns Tuple of 2D-arrays (Hp,GS):
        Hp: Modified H: permutation of columns (The code doesn't change)
        tGS: Transposed Systematic Coding matrix associated to Hp.

    """

    H = np.copy(MATRIX)
    m,n = H.shape
    
    if n>100 and use_sparse:
        sparse = 1
    else:
        sparse = 0 
        
    P1 = np.identity(n,dtype=int)
    
    Hrowreduced = GaussJordan(H)
    
    k = n - sum([a.any() for a in Hrowreduced ])

    ## After this loop, Hrowreduced will have the form H_ss : | I_(n-k)  A |
    permut = np.array(list(range(n)))

    while(True):
        zeros = [i for i in range(min(m,n)) if not Hrowreduced[i,i]]
        if len(zeros)==0:
            break
        indice_colonne_a = min(zeros)
        list_ones = [j for j in range(indice_colonne_a+1,n) if Hrowreduced[indice_colonne_a,j] ]
        if not len(list_ones):
            break

        indice_colonne_b = min(list_ones)
        
        aux = np.copy(Hrowreduced[:,indice_colonne_a])
        Hrowreduced[:,indice_colonne_a] = Hrowreduced[:,indice_colonne_b]
        Hrowreduced[:,indice_colonne_b] = aux 
        
        aux = np.copy(P1[:,indice_colonne_a])
        P1[:,indice_colonne_a] = P1[:,indice_colonne_b]
        P1[:,indice_colonne_b] = aux
        
    ############ NOW, Hrowreduced has the form: | I_(n-k)  A | , the permutation above makes it look like : 
    ########### |A  I_(n-k)|
    
    P1 = P1.T
    identity = list(range(n))
    sigma = identity[n-k:]+identity[:n-k]
    
    P2 = np.zeros(shape=(n,n),dtype=int)
    P2[identity,sigma] = np.ones(n)
    
    if sparse:
        P1 = csr_matrix(P1)
        P2 = csr_matrix(P2)
        H = csr_matrix(H)

    P = BinaryProduct(P2,P1)
    
    if sparse:
        P = csr_matrix(P)
        
    Hp = BinaryProduct(H,np.transpose(P))

    GS = np.zeros((k,n),dtype=int)
    GS[:,:k] = np.identity(k)
    GS[:,k:] = np.transpose(Hrowreduced[:n-k,n-k:])
    
    
    return Hp,np.transpose(GS)

def Coding(tG,v,SNR):
    """
    
    IMPORTANT: if H is large, tG (transposed coding matrix) can be  scipy.sparse.csr_matrix object to speed up calculations.

    Codes a message v with Coding Matrix G, and sends it through a noisy (default)
    channel. 
    
    G's shape is (k,n). 

    Message v is passed to tG: d = tG.tv d is a n-vector turned into a BPSK modulated
    vector x. Then Additive White Gaussian Noise (AWGN) is added. 
    
    SNR is the Signal-Noise Ratio: SNR = 10log(1/variance) in decibels, where variance is the variance of the AWGN.
    Remember: 
    
        1. d = v.G (or (td = tG.tv))
        2. x = BPSK(d) (or if you prefer the math: x = pow(-1,d) )
        3. y = x + AWGN(0,snr) 

    
    Parameters:

    tG: 2D-Array (OR scipy.sparse.csr_matrix)Transposed Coding Matrix obtained from CodingMatrix functions.
    v: 1D-Array, k-vector (binary of course ..) 
    SNR: Signal-Noise-Ratio: SNR = 10log(1/variance) in decibels. 

    -------------------------------

    Returns y

    """
    n,k = tG.shape

    if len(v)!= k:
        raise ValueError(" Size of message v must be equal to number of Coding Matrix G's rows " )  
    
    d = BinaryProduct(tG,v)
    x=pow(-1,d)

    sigma = 10**(-SNR/20)
    e = np.random.normal(0,sigma,size=n)


    y=x+e

    return y

pi = math.pi
def pos(x):
    if  x>0:
        return +1
    elif x<0:
        return -1
    else:
        return 0
def f1(y,sigma):
    """ Normal Density N(1,sigma) """ 
    return(np.exp(-.5*pow((y-1)/sigma,2))/(sigma*math.sqrt(2*pi)))

def fM1(y,sigma):
    """ Normal Density N(-1,sigma) """ 

    return(np.exp(-.5*pow((y+1)/sigma,2))/(sigma*math.sqrt(2*pi)))

def Bits2i(H,i):
    """
    Computes list of elements of N(i)-j:
    List of variables (bits) connected to Parity node i.

    """
    if type(H)!=scipy.sparse.csr_matrix:
        m,n=H.shape
        return ([a for a in range(n) if H[i,a] ])
    
    indj = H.indptr
    indi = H.indices
    
    return [indi[a] for a in range(indj[i],indj[i+1])]


def Nodes2j(tH,j):
    
    """
    Computes list of elements of M(j):
    List of nodes (PC equations) connecting variable j.

    """
    
    return Bits2i(tH,j) 

def InCode(H,x):

    """ Computes Binary Product of H and x. If product is null, x is in the code.

        Returns appartenance boolean. 
    """
        
    return  (BinaryProduct(H,x)==0).all()

def BitsAndNodes(H):
    
    m,n = H.shape
    if type(H)==scipy.sparse.csr_matrix:
        tH = scipy.sparse.csr_matrix(np.transpose(H.toarray()))
    else:
        tH = np.transpose(H)
        
    Bits = [Bits2i(H,i) for i in range(m)]
    Nodes = [Nodes2j(tH,j)for j in range(n)]
    
    return Bits,Nodes
def Decoding_logBP(H,y,SNR,max_iter):
        
    m,n=H.shape

    if not len(y)==n:
        raise ValueError('La taille de y doit correspondre au nombre de colonnes de H')

    if m>=n:
        raise ValueError('H doit avoir plus de colonnes que de lignes')
    
    
    var = 10**(-SNR/10)

    ### ETAPE 0: initialisation 
    
    Lc = 2*y/var
    for k in range(n):
        f.write(str(Lc[k])+'\n')
    Lq=np.zeros(shape=(m,n))
    beta=np.zeros(shape=(d_c*m,1))
    alpha=np.zeros(shape=(d_v*n,1))
    print(alpha)
    Lr = np.zeros(shape=(m,n))
    L = np.zeros(shape=(m,n))
    list=[]
    list2=[]
    list3=[]
    #print(alpha)
    #print(L)
    t=1
    for i in range(m):
        for j in range(n):
            #print(i,j)
            if H[i][j]==1:
                L[i][j]=int(t)
                #print(j)
                list3.extend(str(j))
                t=t+1
    for k in range(0,m*(d_c)):
                    numberlambda.write(str(list3[k])+'\n')
        
    #print(L) 
    count=0
    
    prod=np.prod
    tanh = np.tanh
    log = np.log
    
    #print(H)	
    Bits,Nodes = BitsAndNodes(H)
    while(True):
        l=0
        p=0
        count+=1 #Compteur qui empêche la boucle d'être infinie .. 
        #### ETAPE 1 : Horizontale
        for i in range(m):
            Ni = Bits[i]
            #print(Ni)
            for j in Ni:
                Nij = Ni.copy()

                if j in Nij: Nij.remove(j)
          
                #print(Nij)
                #print(L[i,Nij])
                #print(list)
                list.extend(L[i,Nij]-1)
                #print(list)
                
            
                if count==1:
                   #print(Ni)
                   Lr[i,j] = pos(prod(Lc[Nij]))*min(abs(Lc[Nij]))
                   alpha[l]=Lr[i,j]
                   print(l)
                   
                   l=l+1
                else: 
                    Lr[i,j] = pos(prod(Lq[i,Nij]))*min(abs(Lq[i,Nij]))
                    alpha[l]=Lr[i,j]
                    
                    l=l+1
        print(list)
        print(m)
        print(d_c)
        print(d_v)
        print(n)
        print(L)
        for k in range(0,m*(d_c-1)*(d_c)):
                    numberrow.write(str(list[k])+'\n')
        #### ETAPE 2 : Verticale
        for k in range (1,((m*n)+1)):
            for j in range(n):
                for i in range(m):
                    if(L[i,j]==k):
                        Mj = Nodes[j]
                        #print(Mj)
                        #print("i=",i)
                        Mji = Mj.copy()
                        #print("Mji=",Mji)
                        if i in Mji: Mji.remove(i)
                        #print("Mjinew=",Mji)
                        print("L=",L[Mji,j]-1)
                        list2.extend(L[Mji,j]-1)
                        list2.extend(str(j))
                        Lq[i,j] = Lc[j]+sum(Lr[Mji,j])
                        beta[p]=Lq[i,j]
                        p=p+1
        print (list2)

        #print(list2)
        for k in range(0,(d_v)*(d_v)*n):
            numbercolumn.write(str(list2[k])+'\n')
 
        #### LLR a posterior
        
        L_posteriori = np.zeros(n)
       
        for j in range(n):
            Mj = Nodes[j]

            
            L_posteriori[j] = Lc[j] + sum(Lr[Mj,j])

        x = np.array(L_posteriori <= 0).astype(int)
       
        
            
        product = InCode(H,x)
        if product or count >= max_iter:  
            break
    for k in range(len(alpha)):
        g.write(str(alpha[k]).rstrip('\r\n')+' ')
        h.write(str(beta[k]).rstrip('\r\n')+' ')
    return x

def GaussElimination(MATRIX,B):


    A = np.copy(MATRIX)
    b = np.copy(B)
    n,k = A.shape


    if b.size != n:
        raise ValueError('Size of B must match number of rows of MATRIX to solve MATRIX.X = B')

    for j in range(min(k,n)):
        listeDePivots = [i for i in range(j,n) if A[i,j]]
        if len(listeDePivots)>0:
            pivot = np.min(listeDePivots)
        else:
            continue
        if pivot!=j:
            aux = np.copy(A[j,:])
            A[j,:] = A[pivot,:]
            A[pivot,:] = aux

            aux = np.copy(b[j])
            b[j] = b[pivot]
            b[pivot] = aux

        for i in range(j+1,n):
            if A[i,j]:
                A[i,:] = abs(A[i,:]-A[j,:])
                b[i] = abs(b[i]-b[j])

    return A,b
def DecodedMessage(tG,x):

    n,k = tG.shape 
    
    if n < k:
        raise ValueError('Coding matrix G must have more columns than rows to solve the linear system on v\': G\'v\' = x\'')
    
                         
    rtG, rx = GaussElimination(tG,x)
    
    rank = sum([a.any() for a in rtG])
    if rank!= k:
        raise ValueError('Coding matrix G must have full rank = k to solve G\'v\' = x\'')
            
    message=np.zeros(k).astype(int)

    message[k-1]=rx[k-1]
    for i in reversed(range(k-1)):
        message[i]=abs(rx[i]-BinaryProduct(rtG[i,list(range(i+1,k))],message[list(range(i+1,k))]))

    return message


n = 12  # Number of columns
d_v = 4 # Number of ones per column, must be lower than d_c (because H must have more rows than columns)
d_c = 6 # Number of ones per row, must divide n (because if H has m rows: m*d_c = n*d_v (compute number of ones in H))

H = RegularH(n,d_v,d_c)
#H=np.array([[1,1,1,0,0,0],
 #           [0,0,0,1,1,1],
  #          [1,1,0,0,1,0],
   #         [0,0,1,1,0,1]])


#print("Regular parity-check matrix H({},{},{}):\n\n".format(n,d_v,d_c),H)
tG = CodingMatrix(H)
#print("Transposed Coding Matrix tG that goes with H above is:\n\n",tG)
#print("\n With G,H you can code messages of {} bits into codewords of {} bits because G's shape is {}\n".format(tG.shape[1],tG.shape[0],tG.T.shape))
#new_H,sys_tG = CodingMatrix_systematic(H) #If i'm willing to use only systematic form, I prefer to overwrite
                                                # H and simply compute:
H,tG = CodingMatrix_systematic(H)
#print("The new H is:\n\n",new_H)
#print("\nSystematic tG is:\n\n",sys_tG)

### if you want to check if H.tG = 0:
#print("\nBinary Product HxG' \n\n",BinaryProduct(H,tG))
n,k = tG.shape
#snr = 8
k,n
#v = np.random.randint(2,size=k)
#print (k)
#print("The k-bits message v: \n",v)
#y = Coding(tG,v,snr)
#print (y)
#print("The n-bits message received after transmission:\n\n",y)
#x_decoded = Decoding_logBP(H,y,snr,5)
#print("The decoded n-bits codeword is:\n",x_decoded)
#print("H.x' = ",BinaryProduct(H,x_decoded))
#v_received = DecodedMessage(tG,x_decoded)
#print("The k-bits decoded message v_recieved is:\n",v_received)
#print("\nThe k-bits original message v is:\n",v)
#H=np.array([[1,1,1,0,0,0],
#            [0,0,0,1,1,1],
 #           [1,1,0,0,1,0],
  #          [0,0,1,1,0,1]])
tG = CodingMatrix(H)
#print(tG)
sample_size = 1
input_snr=input(' ')
snrint=int(input_snr)
snr=snrint*0.1
#print(snr)
max_iter = 1
matches=[]
t = time()
for i in range(sample_size): 
    v_sent = np.random.randint(2,size=k)
    print(v_sent)
    for k in range(len(v_sent)):
        o.write(str(v_sent[k])+'\n')
    y = Coding(tG,v_sent,snr)
    for k in range(len(y)):
        p.write(str(y[k])+'\n')
    #y=np.array([0.66379454, -0.59497408, -2.15299592,  2.84935889,  6.53625405,  2.87575803])
    
    x_decoded= Decoding_logBP(H,y,snr,max_iter)
    v_received = DecodedMessage(tG,x_decoded)
    #print(v_received)
    matches.append((v_sent==v_received).all())
t = time()-t
        #print(matches)
        #print("Sucessful decoding rate for snr = {} and max_iter = {} is {}%".format(snr,max_iter,sum(matches)*100/sample_size))
        #print("\nDecoding time of {} messages in seconds:".format(sample_size),t)

