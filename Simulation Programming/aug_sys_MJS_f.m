function dX = aug_sys_MJS_f(t,X,K,L,ifLearned,expl_noise_freq1,expl_noise_freq2,A,B1,B2)
xn=size(A,2);
x = X(1:xn);

if ifLearned==1;   

u = sum(0.1*exp(1)^(-0.01*t)*sin(expl_noise_freq1*t)+0.5*cos(expl_noise_freq1*t),2)-K*x; 
w = sum(0.5*exp(1)^(-0.1*t)*sin(expl_noise_freq2*t)+5*cos(expl_noise_freq2*t),2)+0.1*exp(1)^(-0.01*t)+L*x;


else
    u = -K*x;   
     w = L*x;    
end
dx = A*x+B1*u+B2*w;
dxx = kron(x',x')';
dux = kron(x',u')';
dwx = kron(x',w')';
dX  = [dx;dxx;dux;dwx];
end
