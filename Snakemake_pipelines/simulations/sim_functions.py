import msprime
import numpy as np
import pandas as pd
import random

def missing(inf, outf, miss):

    geno=np.loadtxt(inf,dtype='float', delimiter=",")

    n=int(np.shape(geno)[0] * np.shape(geno)[1] * float(miss))

    index_9 = np.random.choice(geno.size, n, replace=False)
    geno.ravel()[index_9] = 9
    np.savetxt(fname=outf, X=geno, delimiter=",", fmt="%d")



def subset_inds(genof, snpf, indf, subind, s1, s2, s3, s4, e1,
    e2,e3,e4,e5,e6, npop, outg_pc, outg, outi, outs, flag):
    geno=np.loadtxt(genof,dtype='float', delimiter=",")
    snp=np.loadtxt(snpf,dtype='str', delimiter="\t")
    ind=np.loadtxt(indf,dtype='str', delimiter="\t")


    if int(subind) < np.shape(ind)[0]:

        sub1=int(int(subind)/int(npop))
        #gen=np.array([], dtype=np.int64).reshape(len(geno),0)
        sub_p=[]
        for i in range(len(sn)):

            if i<4:
                if flag=='pop':
                    p=np.where(ind[:,2]==sn[i])[0]
                    p1=random.sample(list(p), sub1)
                elif flag=='ind':
                    x1=[x + int(i*len(ind)/int(npop)) for x in list(range(sub1))]
                    p1=x1[:sub1]

            elif i>=4:
                p1=np.where(ind[:,2]==sn[i])[0]

            sub_p.extend(p1)

        sub_p1=np.sort(sub_p)
        geno_out=geno[:,sub_p1]

        #change ind file
        ind_out=ind[sub_p1,:]

        poly=[]
        for x in range(len(geno_out)):
            xxx = geno_out[x,geno_out[x,:]!=9]
            if len(np.unique(xxx)) >=2:
                poly.append(x)

        geno_out1=geno_out[poly,:]

        #make a snp file
        snp1=snp[poly,:]

    elif int(subind) == int(npop):
        geno_out1=geno
        ind_out=ind
        snp1=snp

    np.savetxt(fname=outg_pc, X=(geno_out1).astype(int), delimiter=",", fmt="%d")
    np.savetxt(fname=outg, X=(geno_out1).astype(int), delimiter="", fmt="%s")
    np.savetxt(fname=outi, X=ind_out, delimiter="\t", fmt='%s')
    np.savetxt(fname=outs, X=snp1, delimiter="\t", fmt='%s')


def subset_indsold(genof, snpf, indf, subind, s1, s2, s3, s4, e1,
    e2,e3,e4,e5,e6, npop, outg_pc, outg, outi, outs, flag):
    geno=np.loadtxt(genof,dtype='float', delimiter=",")
    snp=np.loadtxt(snpf,dtype='str', delimiter="\t")
    ind=np.loadtxt(indf,dtype='str', delimiter="\t")

    # subset individusls
    sn=[s1,s2,s3,s4,e1,e2,e3,e4,e5,e6]

    sub1=int(int(subind)/int(npop))
    #gen=np.array([], dtype=np.int64).reshape(len(geno),0)
    sub_p=[]
    for i in range(len(sn)):

        if i<4:
            if flag=='pop':
                p=np.where(ind[:,2]==sn[i])[0]
                p1=random.sample(list(p), sub1)
            elif flag=='ind':
                x1=[x + int(i*len(ind)/int(npop)) for x in list(range(sub1))]
                p1=x1[:sub1]

        elif i>=4:
            p1=np.where(ind[:,2]==sn[i])[0]

        sub_p.extend(p1)

    sub_p1=np.sort(sub_p)
    geno_out=geno[:,sub_p1]

    #change ind file
    ind_out=ind[sub_p1,:]

    poly=[]
    for x in range(len(geno_out)):
        xxx = geno_out[x,geno_out[x,:]!=9]
        if len(np.unique(xxx)) >=2:
            poly.append(x)

    geno_out1=geno_out[poly,:]

    #make a snp file
    snp1=snp[poly,:]


    np.savetxt(fname=outg_pc, X=(geno_out1).astype(int), delimiter=",", fmt="%d")
    np.savetxt(fname=outg, X=(geno_out1).astype(int), delimiter="", fmt="%s")
    np.savetxt(fname=outi, X=ind_out, delimiter="\t", fmt='%s')
    np.savetxt(fname=outs, X=snp1, delimiter="\t", fmt='%s')
