
clc;clear all
%-----------------------------------------------------
A1 = [-0.0500, 6.0000, 0, 0, -0.4775;
      0, -3.4725, 3.4725, 0, 0;
      -5.8775, 0, -13.0210, 0, 0;
      4.0000, 0, 0, 0, 0.7958;
      6.2832, 0, 0, 0, 0];
B1 = [0; 0; 13.0210; 0; 0];
B12 = [-6.0000, 1.0000;
       0, 0;
       0, 0;
       0, -1.6667;
       0, 0];
A2 = [-0.0400, 4.5000, 0, 0, -0.3581;
      0, -3.1565, 3.1565, 0, 0;
      -5.8050, 0, -14.4675, 0, 0;
      4.0000, 0, 0, 0, 0.7958;
      6.2832, 0, 0, 0, 0];
B2 = [0; 0; 14.4675; 0; 0];
B22 = [-4.5000, 1.0000;
       0, 0;
       0, 0;
       0, -2.2222;
       0, 0];
A3 = [-0.0500, 6.0000, 0, 0, -0.4775;
      0, -3.4725, 3.4725, 0, 0;
      -5.8775, 0, -13.0210, 0, 0;
      4.0000, 0, 0, 0, 0.7958;
      6.2832, 0, 0, 0, 0];
B3 = [0; 0; 13.0210; 0; 0];
B32 = [-6.0000, 1.0000;
       0, 0;
       0, 0;
       0, -1.6667;
       0, 0];
A4 = [-0.0400, 4.5000, 0, 0, -0.3581;
      0, -3.1565, 3.1565, 0, 0;
      -5.8050, 0, -14.4675, 0, 0;
      4.0000, 0, 0, 0, 0.7958;
      6.2832, 0, 0, 0, 0];
B4 = [0; 0; 14.4675; 0; 0];
B42 = [-4.5000, 1.0000;
       0, 0;
       0, 0;
       0, -2.2222;
       0, 0];
xn = size(A1,2);un = size(B1,2);wn = size(B12,2);

% [A1, B1, B12] = system_model_revise(1); 
% 
% [A2, B2, B22] = system_model_revise(2); 
% [A3, B3, B32] = system_model_revise(1); 
% [A4, B4, B42] = system_model_revise(2); 

Tr_Q01 =[-1.6 1.6;2.0 -2.0];
Tr_Q02 =[-3.6 3.6;1.0 -1.0];

Tr_Q_0=kron(Tr_Q01, eye(size(Tr_Q02,2))) + kron(eye(size(Tr_Q01,2)), Tr_Q02);
Tr_Q=Tr_Q_0
size(Tr_Q)


A01 = A1+0.5*eye(xn)*Tr_Q(1,1); e1=eig(A01);
A02 = A2+0.5*eye(xn)*Tr_Q(2,2); e2=eig(A02);
A03 = A3+0.5*eye(xn)*Tr_Q(3,3); e3=eig(A03);
A04 = A4+0.5*eye(xn)*Tr_Q(4,4); e4=eig(A04);

e={e1,e2,e3,e4}
for i=1:4
real(e{i})
end

Q01 = 0.1*eye(xn);
R01 = 0.1;  

Q=Q01;
R=R01; 
gamma=3.5;

x010 = [0.1;-0.1;0.1;-0.1;0.1]; 
x0=x010;
size(x0)

load ('multi-Markov_results.mat');
%% 
Initial_K = [0 0 -1.15 0 0];
Initial_L = zeros(wn,xn);

K1 = Initial_K;L1 = Initial_L;
K2 = Initial_K;L2 = Initial_L;
K3= Initial_K;L3 = Initial_L;
K4= Initial_K;L4 = Initial_L;
A01 = A1+0.5*eye(xn)*Tr_Q(1,1);
A02 = A2+0.5*eye(xn)*Tr_Q(2,2);
A03 = A3+0.5*eye(xn)*Tr_Q(3,3);
A04 = A4+0.5*eye(xn)*Tr_Q(4,4);
eig(A01-B1*K1+B12*L1)
eig(A02-B2*K2+B22*L2)
eig(A03-B3*K3+B32*L3)
eig(A04-B4*K4+B42*L4)
%% 
ifLearned=1;   
N=500 ;   
T1  =Jump_time(2)/100;          
T2  =(Jump_time(3)-Jump_time(2))/100;          
T3  =(Jump_time(4)-Jump_time(3))/100;   
T4  =(Jump_time(5)-Jump_time(4))/100;

expl_noise_freq1 = (rand(un,100)-.5)*100; 
expl_noise_freq2 = (rand(wn,100)-.5)*100; 

Dxx1=[];Dxx2=[];Dxx3=[];Dxx4=[];
Ixx1=[];Ixx2=[];Ixx3=[];Ixx4=[];
Ixu1=[];Ixu2=[];Ixu3=[];Ixu4=[];
Ixw1=[];Ixw2=[];Ixw3=[];Ixw4=[];
x01=x0;
X1=[x01;kron(x01',x01')';kron(x01,zeros(un,1));kron(x01,zeros(wn,1))]';
x1_save=[];x2_save=[];t1_save=[];t2_save=[];x3_save=[];t3_save=[];x4_save=[];t4_save=[];

for i=1:N
    [t1,X1] = ode45(@(t,x)aug_sys_MJS_f(t,x,K1,L1,ifLearned,expl_noise_freq1,expl_noise_freq2,A1,B1,B12), [(i-1)*T1,i*T1],X1(end,:)');
    Dxx1=[Dxx1;kron(X1(end,1:xn),X1(end,1:xn))-kron(X1(1,1:xn),X1(1,1:xn))];
    Ixx1=[Ixx1;X1(end,xn+1:xn+xn^2)-X1(1,xn+1:xn+xn^2)];
    Ixu1=[Ixu1;X1(end,xn+xn^2+1:xn+xn^2+xn*un)-X1(1,xn+xn^2+1:xn+xn^2+xn*un)]; 
    Ixw1=[Ixw1;X1(end,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)-X1(1,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)]; 
    x1_save=[x1_save;X1];
    t1_save=[t1_save;t1];
end
x02=X1(end,1:xn)';
X2=[x02;kron(x02',x02')';kron(x02,zeros(un,1));kron(x02,zeros(wn,1))]';

for i=1:N
    [t2,X2] = ode45(@(t,x)aug_sys_MJS_f(t,x,K2,L2,ifLearned,expl_noise_freq1,expl_noise_freq2,A2,B2,B22), [(i-1)*T2,i*T2],X2(end,:)');
    Dxx2=[Dxx2;kron(X2(end,1:xn),X2(end,1:xn))-kron(X2(1,1:xn),X2(1,1:xn))];
    Ixx2=[Ixx2;X2(end,xn+1:xn+xn^2)-X2(1,xn+1:xn+xn^2)];
    Ixu2=[Ixu2;X2(end,xn+xn^2+1:xn+xn^2+xn*un)-X2(1,xn+xn^2+1:xn+xn^2+xn*un)]; 
    Ixw2=[Ixw2;X2(end,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)-X2(1,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)]; 
    x2_save=[x2_save;X2];
    t2_save=[t2_save;t2];
end
x03=X2(end,1:xn)';
X3=[x03;kron(x03',x03')';kron(x03,zeros(un,1));kron(x03,zeros(wn,1))]';
for i=1:N
    [t3,X3] = ode45(@(t,x)aug_sys_MJS_f(t,x,K3,L3,ifLearned,expl_noise_freq1,expl_noise_freq2,A3,B3,B32), [(i-1)*T3,i*T3],X3(end,:)');
    Dxx3=[Dxx3;kron(X3(end,1:xn),X3(end,1:xn))-kron(X3(1,1:xn),X3(1,1:xn))];
    Ixx3=[Ixx3;X3(end,xn+1:xn+xn^2)-X3(1,xn+1:xn+xn^2)];
    Ixu3=[Ixu3;X3(end,xn+xn^2+1:xn+xn^2+xn*un)-X3(1,xn+xn^2+1:xn+xn^2+xn*un)]; 
    Ixw3=[Ixw3;X3(end,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)-X3(1,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)]; 
    x3_save=[x3_save;X3];
    t3_save=[t3_save;t3];
end
x04=X3(end,1:xn)';
X4=[x04;kron(x04',x04')';kron(x04,zeros(un,1));kron(x04,zeros(wn,1))]';
for i=1:N
    [t4,X4] = ode45(@(t,x)aug_sys_MJS_f(t,x,K4,L4,ifLearned,expl_noise_freq1,expl_noise_freq2,A4,B4,B42), [(i-1)*T4,i*T4],X4(end,:)');
    Dxx4=[Dxx4;kron(X4(end,1:xn),X4(end,1:xn))-kron(X4(1,1:xn),X4(1,1:xn))];
    Ixx4=[Ixx4;X4(end,xn+1:xn+xn^2)-X4(1,xn+1:xn+xn^2)];
    Ixu4=[Ixu4;X4(end,xn+xn^2+1:xn+xn^2+xn*un)-X4(1,xn+xn^2+1:xn+xn^2+xn*un)]; 
    Ixw4=[Ixw4;X4(end,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)-X4(1,xn+xn^2+xn*un+1:xn+xn^2+xn*un+xn*wn)]; 
    x4_save=[x4_save;X4];
    t4_save=[t4_save;t4];
end
  save('mode4_samping_data.mat','Dxx1','Ixx1','Ixu1','Ixw1','Dxx2','Ixx2','Ixu2','Ixw2','Dxx3','Ixx3','Ixu3','Ixw3','Dxx4','Ixx4','Ixu4','Ixw4')

load('area1_Initial.mat')
% [Initial_P1,Initial_K1,Initial_L1,Initial_P2,Initial_K2,Initial_L2,Initial_P3,Initial_K3,Initial_L3,Initial_P4,Initial_K4,Initial_L4]=...
%     Diao_zero_sum4(Q,R,Tr_Q,A1,A2,A3,A4,B1,B12,B2,B22,B3,B32,B4,B42,K1,K2,K3,K4,L1,L2,L3,L4,gamma)
eig(Initial_P1)
eig(Initial_P2)
eig(Initial_P3)
eig(Initial_P4)

P1=Initial_P1;
P2=Initial_P2;
P3=Initial_P3;
P4=Initial_P4;

K1=Initial_K1;
K2=Initial_K2;
K3=Initial_K3;
K4=Initial_K4;

L1=Initial_L1;
L2=Initial_L2;
L3=Initial_L3;
L4=Initial_L4;

eig(A01-B1*K1+B12*L1)
eig(A02-B2*K2+B22*L2)
eig(A03-B3*K3+B32*L3)
eig(A04-B4*K4+B42*L4)


[Opt_P1,Opt_K1,Opt_L1,Opt_P2,Opt_K2,Opt_L2,Opt_P3,Opt_K3,Opt_L3,Opt_P4,Opt_K4,Opt_L4]=...
    Diao_Opt_Matrix4(A1,B1,B12,A2,B2,B22,A3,B3,B32,A4,B4,B42,Q,R,Tr_Q,Initial_P1,Initial_P2,Initial_P3,Initial_P4,K1,K2,K3,K4,L1,L2,L3,L4,gamma)
Real_opt_K1=Opt_K1
Real_opt_L1=Opt_L1
Real_opt_P1=Opt_P1

Real_opt_K2=Opt_K2
Real_opt_L2=Opt_L2
Real_opt_P2=Opt_P2

Real_opt_K3=Opt_K3
Real_opt_L3=Opt_L3
Real_opt_P3=Opt_P3

Real_opt_K4=Opt_K4
Real_opt_L4=Opt_L4
Real_opt_P4=Opt_P4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1=Initial_P1;
P2=Initial_P2;
P3=Initial_P3;
P4=Initial_P4;
%% 
P_old1 = eye(xn);
P_old2 = eye(xn);
P_old3 = eye(xn);
P_old4 = eye(xn);
it = 0;            
p1_save = []; p2_save = [];  p3_save = [];   p4_save = [];   
k1_save = []; k2_save = [];  k3_save = [];   k4_save = [];  
l1_save = []; l2_save = [];  l3_save = [];   l4_save = []; 

k1_save  = norm(K1-Real_opt_K1);  
k2_save  = norm(K2-Real_opt_K2);
k3_save  = norm(K3-Real_opt_K3);
k4_save  = norm(K4-Real_opt_K4);

l1_save  = norm(L1-Real_opt_L1);  
l2_save  = norm(L2-Real_opt_L2);
l3_save  = norm(L3-Real_opt_L3);
l4_save  = norm(L4-Real_opt_L4);

MaxIteration = 300;
P1_mat1 =[];P1_mat2 =[];P1_mat3 =[];P1_mat4 =[];
P1_mat5 =[];P1_mat6 =[];P1_mat7 =[];P1_mat8 =[];
P1_mat9 =[];P1_mat10 =[];

P2_mat1 =[];P2_mat2 =[];P2_mat3 =[];P2_mat4 =[];
P2_mat5 =[];P2_mat6 =[];P2_mat7 =[];P2_mat8 =[];
P2_mat9 =[];P2_mat10 =[];

P3_mat1 =[];P3_mat2 =[];P3_mat3 =[];P3_mat4 =[];
P3_mat5 =[];P3_mat6 =[];P3_mat7 =[];P3_mat8 =[];
P3_mat9 =[];P3_mat10 =[];

P4_mat1 =[];P4_mat2 =[];P4_mat3 =[];P4_mat4 =[];
P4_mat5 =[];P4_mat6 =[];P4_mat7 =[];P4_mat8 =[];
P4_mat9 =[];P4_mat10 =[];
P1_mat =[];P2_mat =[];P3_mat =[];P4_mat =[];

norm1=norm(P1 - P_old1, 2);
norm2=norm(P2 - P_old2, 2);
norm3=norm(P3 - P_old3, 2);
norm4=norm(P4 - P_old4, 2);

norm_sum =[norm1,norm2,norm3,norm4];
while   max(norm_sum)>1e-8 & it<MaxIteration
    it = it+1                     
    P_old1 = P1;  
    P_old2 = P2;  
    P_old3 = P3; 
    P_old4 = P4; 

    Theta1 = [Dxx1+1*Tr_Q(1,1)*Ixx1,...
     -2*Ixx1*kron(eye(xn),K1'*R)-2*Ixu1*kron(eye(xn),R),...
     (-2*gamma^2)*(Ixw1-Ixx1*(kron(eye(xn),L1')))]; 
  rank(Theta1);
 Theta2 = [Dxx2+1*Tr_Q(2,2)*Ixx2,...
     -2*Ixx2*kron(eye(xn),K2'*R)-2*Ixu2*kron(eye(xn),R),...
     (-2*gamma^2)*(Ixw2-Ixx2*(kron(eye(xn),L2')))];   
   rank(Theta2);

 Theta3 = [Dxx3+1*Tr_Q(3,3)*Ixx3,...
     -2*Ixx3*kron(eye(xn),K3'*R)-2*Ixu3*kron(eye(xn),R),...
     (-2*gamma^2)*(Ixw3-Ixx3*(kron(eye(xn),L3')))];   
   rank(Theta3);

  Theta4 = [Dxx4+1*Tr_Q(4,4)*Ixx4,...
     -2*Ixx4*kron(eye(xn),K4'*R)-2*Ixu4*kron(eye(xn),R),...
     (-2*gamma^2)*(Ixw4-Ixx4*(kron(eye(xn),L4')))];   
   rank(Theta4);

D1 =Q+Tr_Q(1,2)*P2+Tr_Q(1,3)*P3+Tr_Q(1,4)*P4+K1'*R*K1-gamma^2*L1'*L1;
D2 =Q+Tr_Q(2,1)*P1+Tr_Q(2,3)*P3+Tr_Q(2,4)*P4+K2'*R*K2-gamma^2*L2'*L2;
D3 =Q+Tr_Q(3,1)*P1+Tr_Q(3,2)*P2+Tr_Q(3,4)*P4+K3'*R*K3-gamma^2*L3'*L3;
D4 =Q+Tr_Q(4,1)*P1+Tr_Q(4,2)*P2+Tr_Q(4,3)*P3+K4'*R*K4-gamma^2*L4'*L4;


 Xi1 = -Ixx1*D1(:); 
 Xi2 = -Ixx2*D2(:);      
 Xi3 = -Ixx3*D3(:); 
 Xi4 = -Ixx4*D4(:); 

 pp1 = pinv(Theta1)*Xi1;  
 pp2 = pinv(Theta2)*Xi2; 
 pp3 = pinv(Theta3)*Xi3; 
 pp4 = pinv(Theta4)*Xi4; 

 P1  = reshape(pp1(1:xn*xn), [xn, xn]);
 P2  = reshape(pp2(1:xn*xn), [xn, xn]); 
 P3  = reshape(pp3(1:xn*xn), [xn, xn]); 
 P4  = reshape(pp4(1:xn*xn), [xn, xn]); 

  P1 =0.5*(P1+P1');
  eig(P1);
  P2 =0.5*(P2+P2');
  eig(P2);
  P3 =0.5*(P3+P3');
  eig(P3);
  P4 =0.5*(P4+P4');
  eig(P4);
  

 K1 = reshape(pp1(xn*xn+1:xn*(xn+un)),[un,xn]);
 K2 = reshape(pp2(xn*xn+1:xn*(xn+un)),[un,xn]);
 K3 = reshape(pp3(xn*xn+1:xn*(xn+un)),[un,xn]);
 K4 = reshape(pp4(xn*xn+1:xn*(xn+un)),[un,xn]);

 L1 = reshape(pp1(xn*(xn+un)+1:end),[wn,xn]);
 L2 = reshape(pp2(xn*(xn+un)+1:end),[wn,xn]); 
 L3 = reshape(pp3(xn*(xn+un)+1:end),[wn,xn]); 
 L4 = reshape(pp4(xn*(xn+un)+1:end),[wn,xn]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 p1_save = [p1_save,norm(P1-Real_opt_P1)];   
 k1_save = [k1_save,norm(K1-Real_opt_K1)];    
 p2_save = [p2_save,norm(P2-Real_opt_P2)];  
 k2_save = [k2_save,norm(K2-Real_opt_K2)];    
 l1_save  =[l1_save,norm(L1-Real_opt_L1)];  
 l2_save  =[l2_save,norm(L2-Real_opt_L2)];  
 p3_save = [p3_save,norm(P3-Real_opt_P3)];   
 k3_save = [k3_save,norm(K3-Real_opt_K3)];   
 l3_save  =[l3_save,norm(L3-Real_opt_L3)];  
 p4_save = [p4_save,norm(P4-Real_opt_P4)];   
 k4_save = [k4_save,norm(K4-Real_opt_K4)];   
 l4_save  =[l4_save,norm(L4-Real_opt_L4)];  

norm1=norm(P1-P_old1);
norm2=norm(P2-P_old2);
norm3=norm(P3-P_old3);
norm4=norm(P4-P_old4);
norm_sum =[norm1,norm2,norm3,norm4];
end

P1_mat=[P1_mat1;P1_mat2;P1_mat3;P1_mat4;P1_mat5;P1_mat6;P1_mat7;P1_mat8;P1_mat9;P1_mat10];
P2_mat=[P2_mat1;P2_mat2;P2_mat3;P2_mat4;P2_mat5;P2_mat6;P2_mat7;P2_mat8;P2_mat9;P2_mat10];
P3_mat=[P3_mat1;P3_mat2;P3_mat3;P3_mat4;P3_mat5;P3_mat6;P3_mat7;P3_mat8;P3_mat9;P3_mat10];
P4_mat=[P4_mat1;P4_mat2;P4_mat3;P4_mat4;P4_mat5;P4_mat6;P4_mat7;P4_mat8;P4_mat9;P4_mat10];


a1_p1_save=p1_save;a1_p2_save=p2_save;a1_p3_save=p3_save;a1_p4_save=p4_save;
a1_k1_save=k1_save;a1_k2_save=k2_save;a1_k3_save=k3_save;a1_k4_save=k4_save;
a1_l1_save=l1_save;a1_l2_save=l2_save;a1_l3_save=l3_save;a1_l4_save=l4_save;

save ('area1_pkl_data.mat', 'a1_p1_save', 'a1_p2_save', 'a1_p3_save', 'a1_p4_save','a1_k1_save',...
    'a1_k2_save', 'a1_k3_save', 'a1_k4_save','a1_l1_save', 'a1_l2_save', 'a1_l3_save', 'a1_l4_save')


My_opt_K1=K1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
Real_opt_K1=Real_opt_K1
My_opt_K2=K2
Real_opt_K2=Real_opt_K2
My_opt_K3=K3
Real_opt_K3=Real_opt_K3
My_opt_K4=K4
Real_opt_K4=Real_opt_K4

My_opt_L1=L1
Real_opt_L1=Real_opt_L1
My_opt_L2=L2
Real_opt_L2=Real_opt_L2
My_opt_L3=L3
Real_opt_L3=Real_opt_L3
My_opt_L4=L4
Real_opt_L4=Real_opt_L4

My_opt_P1=P1
Real_opt_P1=Real_opt_P1
My_opt_P2=P2
Real_opt_P2=Real_opt_P2
My_opt_P3=P3
Real_opt_P3=Real_opt_P3
My_opt_P4=P4
Real_opt_P4=Real_opt_P4

norm01=norm(My_opt_P1-Real_opt_P1,2)
norm02=norm(My_opt_P2-Real_opt_P2,2)
norm03=norm(My_opt_P3-Real_opt_P3,2)
norm04=norm(My_opt_P4-Real_opt_P4,2)
%
save('area1_mode4_opt_data.mat','My_opt_P1','Real_opt_P1','My_opt_P2','Real_opt_P2',...
    'My_opt_P3','Real_opt_P3','My_opt_P4','Real_opt_P4','norm01','norm02','norm03','norm04')












