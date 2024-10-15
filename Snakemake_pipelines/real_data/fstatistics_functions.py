import numpy as np
import pandas as pd

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



def sd_ppca(ppcaf, outfile_mu, outfile_std, scale, filter1f):
    filter1=np.loadtxt(filter1f,dtype='int', delimiter=",")

    f3_0=np.loadtxt(ppcaf[0],dtype='float', delimiter=",").reshape(1,-1)
    for i in range(1,len(ppcaf)):
        f3_0=np.concatenate((f3_0, np.loadtxt(ppcaf[i],dtype='float', delimiter=",").reshape(1,-1)),axis=0)

    mu=np.mean(f3_0,0)
    var = var = np.array([np.sum( (np.sum(filter1) - filter1) / filter1 * (f3_0[:,k]-mu[k])**2) / len(filter1) for k in range(len(mu))])
    std = var**0.5



    np.savetxt(fname=outfile_mu, X=mu, delimiter=",", fmt="%f")
    np.savetxt(fname=outfile_std, X=std, delimiter=",", fmt="%f")


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

        f2_res = np.sum(f2(datapop, popin[0], popin[1]))
        f3_res = np.sum(f3(datapop, popin[1], popin[0], popin[2]))
        f4_res = np.sum(f4(datapop, popin[0], popin[1], popin[2], popin[3]))

        f2mat=np.zeros([len(popin), len(popin)])
        for i in range(len(popin)-1):
            for j in range(i+1, len(popin)):
                f2mat[i,j] = np.sum(f2(datapop, popin[i], popin[j]))
                f2mat[j,i] = f2mat[i,j]

    elif flag=="pca_methods":

        f2_res = [norm_scale * 4 * np.sum(f2(datapop[:x], popin[0], popin[1])) for x in list(range(1,np.shape(datapop)[0]+1))]
        f3_res = [norm_scale * 4 * np.sum(f3(datapop[:x], popin[1], popin[0], popin[2])) for x in list(range(1,np.shape(datapop)[0]+1))]
        f4_res = [norm_scale * 4 * np.sum(f4(datapop[:x], popin[0], popin[1], popin[2], popin[3])) for x in list(range(1,np.shape(datapop)[0]+1))]

        f2mat=np.zeros([len(popin), len(popin)])
        for i in range(len(popin)-1):
            for j in range(i+1, len(popin)):
                f2mat[i,j] = [norm_scale * 4 * np.sum(f2(datapop[:x], popin[i], popin[j])) for x in list(range(1,np.shape(datapop)[0]+1))][int(scale)]
                f2mat[j,i] = f2mat[i,j]



    elif flag=="adjusted":
        f2mat = adjustedf2(pmat=datapop, popindex=popin, narray=narray)
        f2_res = f2mat[popin[0], popin[1]]
        f3_res = (f2mat[popin[0], popin[1]] + f2mat[popin[1], popin[2]] - f2mat[popin[0], popin[2]])/2
        f4_res = (f2mat[popin[0], popin[3]] + f2mat[popin[1], popin[2]] - f2mat[popin[0], popin[2]] - f2mat[popin[1], popin[3]])/2




    f_res=np.array([f2_res, f3_res, f4_res])

    np.savetxt(fname=outfile, X=f_res, delimiter=",")
    np.savetxt(fname=f2out, X=f2mat, delimiter=",")



def f3stats(dataf, genof, popf, filter1, f3out, scale, f,f2altai_vin):
    f=int(f)
    geno=np.loadtxt(genof,dtype='int', delimiter=",")
    fil1=np.loadtxt(filter1,dtype='int', delimiter=",")
    data=np.loadtxt(dataf,dtype='float', delimiter=",")[:2,:]
    pop=np.loadtxt(popf,dtype='str', delimiter="\t")
    norm_scale=len(geno) - fil1[f-1]

    pop1=np.where(pop[:,0]=='Altai')[0][0]
    pop3 = np.where(pop[:,0]=='Vindija33.19')[0][0]

    pop2_1 = np.where(pop[:,0]=='Les_Cottes_L35MQ25')[0][0]
    pop2_2 = np.where(pop[:,0]=='Goyet_L35MQ25')[0][0]
    pop2_3 = np.where(pop[:,0]=='Mezmaiskaya1_L35MQ25')[0][0]
    pop2_4 = np.where(pop[:,0]=='Mezmaiskaya2_L35MQ25')[0][0]
    pop2_5 = np.where(pop[:,0]=='VindijaG1_L35MQ25')[0][0]
    pop2_6 = np.where(pop[:,0]=='Spy_L35MQ25')[0][0]

    f3_res = [ 1 * np.sum(f3(data, pop1, pop2_1, pop3)),  1 * np.sum(f3(data, pop1, pop2_2, pop3)), 1 * np.sum(f3(data, pop1, pop2_3, pop3)), 1 * np.sum(f3(data, pop1, pop2_4, pop3)), 1 * np.sum(f3(data, pop1, pop2_5, pop3)), 1 * np.sum(f3(data, pop1, pop2_6, pop3))]

    f2val=np.sum(f2(data,pop1,pop3))
    with open(f2altai_vin, 'w') as file:
        file.write(str(f2val))


    np.savetxt(fname=f3out, X=f3_res, delimiter=",")


def f2altvin(f2,f2out):
    #with open(f2[0], 'r') as file:
        #xx0 = float(file.read())
    xx0=[]

    for i in range(len(f2)):
        with open(f2[i], 'r') as file:
            xx1 = float(file.read())
            xx0.append(xx1)

    print(np.mean(xx0))
    print(np.std(xx0))

    np.savetxt(fname=f2out, X=xx0, delimiter=",")



def sd_ppca_all(ppcaf, outfile_mu, outfile_std, scale, filter1f, npcs):
    npcs=int(npcs)
    filter1=np.loadtxt(filter1f,dtype='int', delimiter=",")

    ppca_0=np.loadtxt(ppcaf[0],dtype='float', delimiter=",")[:npcs,:]
    sig=np.sign(ppca_0)

    for i in range(1,len(ppcaf)):
        xx=np.loadtxt(ppcaf[i],dtype='float', delimiter=",")[:npcs,:]

        for b in range(npcs):
            if np.sum(np.sign(xx[b,:])!=sig[b]) / len(sig[b]) >0.5 :
                xx[b,:] = (-1) * xx[b,:]

        ppca_0=np.dstack((ppca_0, xx))

    mu=np.mean(ppca_0,2)

    var=np.zeros(np.shape(mu))
    for j in range(np.shape(mu)[0]):
        for k in range(np.shape(mu)[1]):
            inst1=ppca_0[j,k,:]
            var[j,k] = np.sum( (np.sum(filter1) - filter1) / filter1 * (inst1-mu[j,k])**2) / len(filter1)

    std = var**0.5


    np.savetxt(fname=outfile_mu, X=mu, delimiter=",", fmt="%f")
    np.savetxt(fname=outfile_std, X=std, delimiter=",", fmt="%f")



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


def f2stats(dataf, popf, filter1f, scale, fout_std, fout_mean):

    scale=int(scale)
    pop=np.loadtxt(popf,dtype='str', delimiter="\t")[:,0]
    filter1=np.loadtxt(filter1f,dtype='int', delimiter=",")

    df_array=[]

    for f in range(len(dataf)):

        data=np.loadtxt(dataf[f],dtype='float', delimiter=",")[:scale,:]

        df = pd.DataFrame(index=pop, columns=pop)
        #df_std = pd.DataFrame(index=pop[:,0], columns=pop[:,0])

        for i in range(len(pop)):
            for j in range(len(pop)):

                if i==j:
                    df.loc[pop[i],pop[j]] = 0
                else:
                    df.loc[pop[i],pop[j]] = np.sum(f2(data,i,j))

        df_array.append(df)

    stacked_data = np.stack([df_.values for df_ in df_array], axis=0)
    mean_data = np.mean(stacked_data, axis=0)
    mean_df = pd.DataFrame(mean_data, index=df.index, columns=df.columns)

    var=np.zeros(np.shape(mean_data))
    for j in range(np.shape(mean_data)[0]):
        for k in range(np.shape(mean_data)[1]):
            inst1=stacked_data[:,j,k]
            var[j,k] = np.sum( (np.sum(filter1) - filter1) / filter1 * (inst1-mean_data[j,k])**2) / len(filter1)

    std = var**0.5
    std_df = pd.DataFrame(std, index=df.index, columns=df.columns)

    with pd.option_context('display.max_rows', len(mean_df.index), 'display.max_columns', len(mean_df.columns)):
        mean_df.to_csv(fout_mean, sep=',')

    with pd.option_context('display.max_rows', len(std_df.index), 'display.max_columns', len(std_df.columns)):
        std_df.to_csv(fout_std, sep=',')
