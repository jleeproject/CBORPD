function F=Ttrans(x,Eg)
F=sign(x).*double((abs(x)-Eg)>0).*(abs(x)-Eg);
end