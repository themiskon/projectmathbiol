syms eps Cmax bz Ks P Z dz g H N Mp K 
eq1=eps*Cmax*(bz*P/(K+bz*P)*Z)-dz*Z;
eq2=(g*N/(H+N)-Mp)*P-Cmax*(bz*P/(K+bz*P))*Z;
eq3=-(g*N/(H+N)-Mp)*P;

sol= solve([eq1==0,eq2==0,eq3==0],[Z, P, N]);
sol.Z
sol.P
sol.N