function [Opt_P1,Opt_K1,Opt_L1,Opt_P2,Opt_K2,Opt_L2,Opt_P3,Opt_K3,Opt_L3,Opt_P4,Opt_K4,Opt_L4]=...
    Diao_Opt_Matrix(A1,B1,B12,A2,B2,B22,A3,B3,B32,A4,B4,B42,Q,R,Tr_Q,Initial_P1,Initial_P2,Initial_P3,Initial_P4,K1,K2,K3,K4,L1,L2,L3,L4,gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xn=size(Q,2);un=size(R,2);
    A01=A1+0.5*Tr_Q(1,1)*eye(xn);
    A02=A2+0.5*Tr_Q(2,2)*eye(xn);
    A03=A3+0.5*Tr_Q(3,3)*eye(xn);
    A04=A4+0.5*Tr_Q(4,4)*eye(xn);
P_old1 = eye(xn);P_old2=eye(xn);P_old3 = eye(xn);P_old4 = eye(xn);

P11 = Initial_P1;P21 = Initial_P2;P31 = Initial_P3;P41 = Initial_P4;


norm1=norm(P11 - P_old1);
norm2=norm(P21 - P_old2);
norm3=norm(P31 - P_old3);
norm4=norm(P41 - P_old4);
norm_sum =[norm1,norm2,norm3,norm4];



j=0;

while  max(norm_sum)>1e-7 
    j=j+1
P_old1 = P11; P_old2 = P21; P_old3 = P31; P_old4 = P41;

B01 = [B1,B12];
m01 = size(B1,2);
m02 = size(B12,2);

R01 = [R zeros(m01,m02);  zeros(m02,m01) -gamma^(2)*eye(m02)];
QK01=Tr_Q(1,2)*P21+Tr_Q(1,3)*P31+Tr_Q(1,4)*P41;
Q01 =Q+QK01;
P11 = care(A01,B01,Q01,R01);

B02 = [B2,B22];
m001 = size(B2,2);
m002 = size(B22,2);

R02 = [R zeros(m001,m002);  zeros(m002,m001) -gamma^(2)*eye(m002)];
QK02=Tr_Q(2,1)*P11+Tr_Q(2,3)*P31+Tr_Q(2,4)*P41;
Q02 =Q+QK02;
P21 = care(A02,B02,Q02,R02);

B03 = [B3,B32];
m031 = size(B3,2);
m032 = size(B32,2);

R03 = [R zeros(m031,m032);  zeros(m032,m031) -gamma^(2)*eye(m032)];
QK03=Tr_Q(3,1)*P11+Tr_Q(3,2)*P21+Tr_Q(3,4)*P41;
Q03 =Q+QK03;
P31 = care(A03,B03,Q03,R03);

B04 = [B4,B42];
m041 = size(B4,2);
m042 = size(B42,2);

R04 = [R zeros(m041,m042);  zeros(m042,m041) -gamma^(2)*eye(m042)];
QK04=Tr_Q(4,1)*P11+Tr_Q(4,2)*P21+Tr_Q(4,3)*P31;
Q04 =Q+QK04;
P41 = care(A04,B04,Q04,R04);


norm1=norm(P11-P_old1);
norm2=norm(P21-P_old2);
norm3=norm(P31-P_old3);
norm4=norm(P41-P_old4);
norm_sum =[norm1,norm2,norm3,norm4];

end
Opt_P1=P11;
Opt_K1=inv(R)*(B1')*Opt_P1;
Opt_L1=gamma^(-2)*(B12')*Opt_P1;

Opt_P2=P21;
Opt_K2=inv(R)*(B2')*Opt_P2;
Opt_L2=gamma^(-2)*(B22')*Opt_P2;

Opt_P3=P31;
Opt_K3=inv(R)*(B3')*Opt_P3;
Opt_L3=gamma^(-2)*(B32')*Opt_P3;

Opt_P4=P41;
Opt_K4=inv(R)*(B4')*Opt_P4;
Opt_L4=gamma^(-2)*(B42')*Opt_P4;
end




