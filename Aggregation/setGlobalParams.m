function [ ] = setGlobalParams( rval,sigval,alphaval,Nval,kval,bval,betaval,mval,Kmaxval )
    global r sig alpha N k b beta m Kmax
    r=rval;
    sig=sigval;
    alpha=alphaval;
    N=Nval;
    k=kval;
    b=bval;
    beta=betaval;
    m=mval;
    Kmax=Kmaxval;
end
