function [ ] = setGlobalParams( rval,sigval,alphaval,Nval,kval,bval,betaval,mval,Kmaxval,sval,indexval,treatmentval,scheduleval )
    global r sig alpha N k b beta m Kmax s index treatment schedule
    r=rval;
    sig=sigval;
    alpha=alphaval;
    N=Nval;
    k=kval;
    b=bval;
    beta=betaval;
    m=mval;
    Kmax=Kmaxval;
    s=sval;
    index=indexval;
    treatment=treatmentval;
    schedule=scheduleval;
end
