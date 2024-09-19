import numpy as np
import pandas as pd


def maketable(f2_1, f2_2, f2_3, f3_1, f3_2, f3_3, f4_1, f4_2, f4_3, truef, outfile):

    tval=np.loadtxt(truef,dtype='float', delimiter=",")/118000
    t12=tval[0,1]
    t34=tval[2,3]
    t23=tval[1,2]
    t14=tval[0,3]
    t13=tval[0,2]
    t24=tval[1,3]

    df2_1=pd.read_csv(f2_1, sep=",", header=0, index_col=False, low_memory=False)
    df2_1.columns=['Method', 'Statistic', 'Mean', 'Standard Error']
    df2_1["Simulation"]="Population-based"
    #df2_1['Statistic']=df2_1['Statistic'].astype(str)

    df2_2=pd.read_csv(f2_2, sep=",", header=0, index_col=False, low_memory=False)
    df2_2.columns=['Method', 'Statistic', 'Mean', 'Standard Error']
    df2_2["Simulation"]="Individual-based"
    #df2_2['Statistic']=df2_2['Statistic'].astype(str)

    df2_3=pd.read_csv(f2_3, sep=",", header=0, index_col=False, low_memory=False)
    df2_3.columns=['Method', 'Statistic', 'Mean', 'Standard Error']
    df2_3["Simulation"]="Individual-based with missing values"
    #df2_3['Statistic']=df2_3['Statistic'].astype(str)

    df3_1=pd.read_csv(f3_1, sep=",", header=0, index_col=False, low_memory=False)
    df3_1.columns=['Method', 'Standard Error', 'Mean', 'Statistic']
    df3_1["Simulation"]="Population-based"

    df3_2=pd.read_csv(f3_2, sep=",", header=0, index_col=False, low_memory=False)
    df3_2.columns=['Method', 'Standard Error', 'Mean', 'Statistic']
    df3_2["Simulation"]="Individual-based"

    df3_3=pd.read_csv(f3_3, sep=",", header=0, index_col=False, low_memory=False)
    df3_3.columns=['Method', 'Standard Error', 'Mean', 'Statistic']
    df3_3["Simulation"]="Individual-based with missing values"

    df4_1=pd.read_csv(f4_1, sep=",", header=0, index_col=False, low_memory=False)
    df4_1.columns=['Method', 'Standard Error', 'Mean', 'Statistic']
    df4_1["Simulation"]="Population-based"

    df4_2=pd.read_csv(f4_2, sep=",", header=0, index_col=False, low_memory=False)
    df4_2.columns=['Method', 'Standard Error', 'Mean', 'Statistic']
    df4_2["Simulation"]="Individual-based"

    df4_3=pd.read_csv(f4_3, sep=",", header=0, index_col=False, low_memory=False)
    df4_3.columns=['Method', 'Standard Error', 'Mean', 'Statistic']
    df4_3["Simulation"]="Individual-based with missing values"


    df=pd.concat([df2_1, df2_2, df2_3, df3_1, df3_2, df3_3, df4_1, df4_2, df4_3],axis=0,join='inner')
    df=df.reset_index(drop=True)


    df['Standard Error'] = df['Standard Error']/118000
    df['Mean'] = df['Mean']/118000


    df.loc[df["Method"] ==  "PCA1_val_scale8","Method"] = "LSE_scale8"
    df.loc[df["Method"] ==  "PCA1_val_scale12","Method"] = "LSE_scale12"
    df.loc[df["Method"] ==  "ppca_direct_val_scale8","Method"] = "PPCA_scale8"
    df.loc[df["Method"] ==  "ppca_direct_val_scale12","Method"] = "PPCA_scale12"
    df.loc[df["Method"] ==  "admixtools2Norm_scale8","Method"] = "ADMIXTOOLS 2"
    df.loc[df["Method"] ==  "admixtools2Norm_scale12","Method"] = "ADMIXTOOLS 2"
    df.loc[df["Method"] ==  "emu_val_scale8","Method"] = "PCA_scale8"
    df.loc[df["Method"] ==  "emu_val_scale12","Method"] = "PCA_scale12"
    df.loc[df["Method"] ==  "ppca_miss_val_scale8","Method"] = "PPCA_scale8"
    df.loc[df["Method"] ==  "ppca_miss_val_scale12","Method"] = "PPCA_scale12"

    df.loc[df["Statistic"] ==  "1;3,4","Statistic"] = "f3(X1;X3,X4)"
    df.loc[df["Statistic"] ==  "4;1,2","Statistic"] = "f3(X4;X1,X2)"
    df.loc[df["Statistic"] ==  "1,3;2,4","Statistic"] = "f4(X1,X3;X2,X4)"


    df.loc[df["Statistic"]==1,"Statistic"] = "f2(X1,X2)"
    df.loc[df["Statistic"]==2,"Statistic"] = "f2(X3,X4)"
    df.loc[df["Statistic"]==3,"Statistic"] = "f2(X2,X3)"
    df.loc[df["Statistic"]==4,"Statistic"] = "f2(X1,X4)"

    df['Statistic'] = df['Statistic'].str.replace("'", "")
    df['Statistic'] = df['Statistic'].str.replace('"', '')
    df['Statistic'] = df['Statistic'].replace('"', '', regex=True)

    df.loc[df['Statistic']=="f2(X1,X2)","True value"] = t12
    df.loc[df['Statistic']=="f2(X2,X3)","True value"] = t23
    df.loc[df['Statistic']=="f2(X3,X4)","True value"] = t34
    df.loc[df['Statistic']=="f2(X1,X4)","True value"] = t14

    df.loc[df['Statistic']=="f3(X1;X3,X4)","True value"] = (t13+t14-t34)/2
    df.loc[df['Statistic']=="f3(X4;X1,X2)","True value"] = (t14+t24-t12)/2

    df.loc[df['Statistic']=="f4(X1,X3;X2,X4)","True value"] = (t13+t24-t12-t34)/2

    df["Bias"] = np.round(df["Mean"] - df["True value"],4)


    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
                df.to_csv(outfile, sep=',', index=False)
