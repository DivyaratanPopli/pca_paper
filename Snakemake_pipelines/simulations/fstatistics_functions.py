import numpy as np


def f2(data, S1, S2):
    return np.square(data[:,S1] - data[:,S2])

def f3(data, S1, S2, S3):
    return (data[:,S1] - data[:,S2]) * (data[:,S1] - data[:,S3])

def f4(data, S1, S2, S3, S4):
    return (data[:,S1] - data[:,S2]) * (data[:,S3] - data[:,S4])
    #return (f2(data, S1, S4) +  f2(data, S2, S3) - f2(data, S1, S2) - f2(data, S3, S4))/2

def truef2(mutation_rate, Ne, st, n, length):
    theta = 2*mutation_rate*length
    ET11 = 2*Ne*(1-(1/n))
    ET22 = 2*Ne*(1-(1/n))
    ET12= 2*Ne*(1-(1/n)) + st
    f2_res = theta * (ET12 - (ET11+ET22)/2)
    return f2_res


def adjustedf2(pmat, popindex, narray):

    f2mat=np.zeros([len(popindex), len(popindex)])
    for j in range(len(popindex)-1):
        for k in range(j+1,len(popindex)):
            pj=pmat[:,j]
            pk=pmat[:,k]
            f2mat[j,k]= np.sum(np.square(pj - pk) - (pj*(1-pj)/(narray[j]-1)) - (pk*(1-pk)/(narray[k]-1)))
            f2mat[k,j]=f2mat[j,k]

    return f2mat


def f4_test(fpos, fppca, popf, S1, S2, S3, S4, scale, fout):
    #ncol=len(np.loadtxt(fpos,dtype='float', delimiter=","))
    #print(fpos, fppca, popf, S1, S2, S3, S4, scale, fout)
    with open(fpos,"r") as f:
            ncol=float(f.read())

    datapop,narray,popin=f_common(dataf=fppca, popf=popf, S1=S1, S2=S2, S3=S3, S4=S4)
    datapop=datapop*2# since there is division by 2 in data2pop function
    xxx=datapop[:int(scale),:].T
    Y=xxx @ xxx.T
    X=np.array([2,2])
    print("Y",Y)
    X=np.array([[Y[0,0]+Y[1,1]-2*Y[0,1] , Y[0,2]+Y[1,3]-Y[0,3]-Y[1,2]] , [Y[0,2]+Y[1,3]-Y[0,3]-Y[1,2] , Y[2,2]+Y[3,3]-2*Y[2,3]]])
    print("X",X)

    #ll=ncol*np.log10((X[0,0]*X[1,1]+X[0,1]**2)/(X[0,0]*X[1,1]))
    ll=ncol*np.log10((X[0,0]*X[1,1])/((X[0,0]*X[1,1]) - X[0,1]**2))

    with open(fout, 'w') as f:
        print(ll,file=f)

def all_ftest(ftest, fadmix, mu, runs, fout):
    flist=[]
    mlist=[]
    rlist=[]
    alist=[]
    i=0
    for m in mu:
        for r in runs:

            val = np.loadtxt(ftest[i],dtype='float', delimiter=",")

            adm= np.loadtxt(fadmix[i],dtype='float', delimiter=",")
            a=adm[2]/adm[3]

            flist.append(val)
            mlist.append(m)
            rlist.append(r)
            alist.append(a)
            i=i+1
    allf=np.stack((mlist,rlist,flist,alist),axis=1)
    np.savetxt(fname=fout, X=allf, delimiter=",")


def ftrue(f2list, flist, outf2, outf):

    f20=np.loadtxt(f2list[0],dtype='float', delimiter=",")
    f0=np.loadtxt(flist[0],dtype='float', delimiter=",")

    for i in range(1,len(f2list)):
        f20=f20+np.loadtxt(f2list[i],dtype='float', delimiter=",")
        f0=f0+np.loadtxt(flist[i],dtype='float', delimiter=",")

    f20=f20/(i+1)
    f0=f0/(i+1)

    np.savetxt(fname=outf2, X=f20, delimiter=",")
    np.savetxt(fname=outf, X=f0, delimiter=",")

def admix_scale(dataf, genof, outfile, f2inf, f2normf):
    data=np.loadtxt(dataf,dtype='float', delimiter=",")[:3] #fourth value is st. dev.
    geno=np.loadtxt(genof,dtype='float', delimiter=",")
    norm_scale=len(geno)
    print(norm_scale, data)

    [f2_res, f3_res, f4_res]=norm_scale * data
    f_res = np.array([f2_res, f3_res, f4_res])

    f2in=np.loadtxt(f2inf,dtype='float', delimiter=",")
    f2norm = f2in * norm_scale

    np.savetxt(fname=outfile, X=f_res, delimiter=",")
    np.savetxt(fname=f2normf, X=f2norm, delimiter=",")

def data2pop(data,sn):

    datapop=np.zeros([len(data), len(sn)])
    narray=[]
    for p in range(len(sn)):
        p_pop=data[:,sn[p]].reshape([len(data),-1])
        narray.append(2*np.shape(p_pop)[1])
        datapop[:,p] = np.nansum(p_pop,1)/(2*np.sum(~np.isnan(p_pop),1))

    popin=list(range(len(sn)))

    return datapop,narray,popin

def f_common(dataf, popf, S1, S2, S3, S4):
    data=np.loadtxt(dataf,dtype='float', delimiter=",")
    pops=np.loadtxt(popf,dtype='str', delimiter="\t")
    data=data.astype(float)
    data[data==9]=np.nan

    s1=np.where(pops[:,2]==S1)[0]
    s2=np.where(pops[:,2]==S2)[0]
    s3=np.where(pops[:,2]==S3)[0]
    s4=np.where(pops[:,2]==S4)[0]
    sn=[s1,s2,s3,s4]

    datapop,narray,popin=data2pop(data=data,sn=sn)
    return datapop,narray,popin



def fstats(dataf, genof, popf, S1, S2, S3, S4, flag, f2out, scale, outfile):

    datapop,narray,popin=f_common(dataf, popf, S1, S2, S3, S4)
    geno=np.loadtxt(genof,dtype='float', delimiter=",")
    norm_scale=len(geno)
    if flag=="noisy":

        f2_res = np.nansum(f2(datapop, popin[0], popin[1]))
        f3_res = np.nansum(f3(datapop, popin[1], popin[0], popin[2]))
        f4_res = np.nansum(f4(datapop, popin[0], popin[1], popin[2], popin[3]))
        print(f2(datapop, popin[0], popin[1]))
        f2mat=np.zeros([len(popin), len(popin)])
        for i in range(len(popin)-1):
            for j in range(i+1, len(popin)):
                f2mat[i,j] = np.nansum(f2(datapop, popin[i], popin[j]))
                f2mat[j,i] = f2mat[i,j]

    elif flag=="pca_methods":

        f2_res = [norm_scale * 4 * np.sum(f2(datapop[:x], popin[0], popin[1])) for x in list(range(1,np.shape(datapop)[0]+1))]
        f3_res = [norm_scale * 4 * np.sum(f3(datapop[:x], popin[1], popin[0], popin[2])) for x in list(range(1,np.shape(datapop)[0]+1))]
        f4_res = [norm_scale * 4 * np.sum(f4(datapop[:x], popin[0], popin[1], popin[2], popin[3])) for x in list(range(1,np.shape(datapop)[0]+1))]

        f2mat=np.zeros([len(popin), len(popin)])
        for i in range(len(popin)-1):
            for j in range(i+1, len(popin)):
                f2mat[i,j] = [norm_scale * 4 * np.sum(f2(datapop[:x], popin[i], popin[j])) for x in list(range(1,np.shape(datapop)[0]+1))][int(scale)-1]
                f2mat[j,i] = f2mat[i,j]



    elif flag=="adjusted":
        f2mat = adjustedf2(pmat=datapop, popindex=popin, narray=narray)
        f2_res = f2mat[popin[0], popin[1]]
        f3_res = (f2mat[popin[0], popin[1]] + f2mat[popin[1], popin[2]] - f2mat[popin[0], popin[2]])/2
        f4_res = (f2mat[popin[0], popin[3]] + f2mat[popin[1], popin[2]] - f2mat[popin[0], popin[2]] - f2mat[popin[1], popin[3]])/2




    f_res=np.array([f2_res, f3_res, f4_res])

    np.savetxt(fname=outfile, X=f_res, delimiter=",")
    np.savetxt(fname=f2out, X=f2mat, delimiter=",")




def make_table(ftrue, fmat, outfile):
    true=np.loadtxt(ftrue,dtype='float', delimiter=",")
    mat=np.loadtxt(fmat,dtype='float', delimiter=",")

    np.fill_diagonal(true,1)
    np.fill_diagonal(mat,1)
    outmat=(mat-true)#*100/true

    np.savetxt(fname=outfile, X=outmat, delimiter=",")

def make_tables(allf, meanf, stdf):
    x0=[]
    for i in range(len(allf)):
        x0.append(np.loadtxt(allf[i], dtype='float', delimiter=","))

    xmean=np.mean(x0,0)
    xstd=np.std(x0,0)

    np.savetxt(fname=meanf, X=xmean, delimiter=",")
    np.savetxt(fname=stdf, X=xstd, delimiter=",")


def f3_f4(allf, meanf, stdf):
    f3_134_all=[]
    f3_412_all=[]
    f4_1324_all=[]
    for i in range(len(allf)):
        f2 = np.loadtxt(allf[i], dtype='float', delimiter=",")
        f3_134 = (f2[0,2] + f2[0,3] - f2[2,3]) /2
        f3_412 = (f2[0,3] + f2[1,3] - f2[0,1]) /2
        f4_1324 = (f2[0,3] + f2[1,2] - f2[0,1] - f2[2,3] ) /2

        f3_134_all.append(f3_134)
        f3_412_all.append(f3_412)
        f4_1324_all.append(f4_1324)

    fmean=[np.mean(f3_134_all), np.mean(f3_412_all), np.mean(f4_1324_all)]
    fstd=[np.std(f3_134_all), np.std(f3_412_all), np.std(f4_1324_all)]

    np.savetxt(fname=meanf, X=fmean, delimiter=",")
    np.savetxt(fname=stdf, X=fstd, delimiter=",")

def f3_f4_true(allf, meanf):


    f2 = np.loadtxt(allf, dtype='float', delimiter=",")
    f3_134_all = (f2[0,2] + f2[0,3] - f2[2,3]) /2
    f3_412_all = (f2[0,3] + f2[1,3] - f2[0,1]) /2
    f4_1324_all = (f2[0,3] + f2[1,2] - f2[0,1] - f2[2,3] ) /2
    fmean=[f3_134_all, f3_412_all, f4_1324_all]


    np.savetxt(fname=meanf, X=fmean, delimiter=",")



def ind2pop(indfile, popfile):

    ind = np.loadtxt(indfile, dtype='str', delimiter="\t")

    ind[:,2]=[ind[i,2].split("_")[0] for i in range(len(ind))]

    np.savetxt(fname=popfile, X=ind, delimiter="\t", fmt="%s")
