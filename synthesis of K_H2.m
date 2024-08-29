m1= 0.128; % mass of the fine stage
m2= 1.977; % mass of the coarse stage
b1= 17.7; % counter electromotive force of VCM
b2= 19.6; % counter electromotive force of linear motor
k1= 6797; % stiffness of the flexure 用测力计测一下
kf1= 17.7;% VCM force constant N/A 
kv1= 0.2; % VCM driver gain A/V
kf2= 24;  % linear PM motor force constant N/A
kv2= 1; % linear PM motor driver gain A/V
f= 4.8; %(unsure)
A=[0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 -k1/m1 -b1/m1 0 k1/m1 b1/m1;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   0 k1/m2 b1/m2 0 -k1/m2 -b2/m2-b1/m2];%之前-b2/m2-b1/m2推错了
B2=[0 0;
    0 0;
    -kf1*kv1/m1 0;
    0 0;
    0 0;
    kf1*kv1/m2 -kf2*kv2/m2];
B1=[0;0;0;0;0;f/m2];%
%  B1=eye(6);
%  B1=eye(6);
C=[100000 0 0 0 0 0;
   0 0 0 1000 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0];
D=[0 0 ;
   0 0 ;
   1 0 ;
   0 1;];
%%%%system in the parameter space
F = [A          -B2;
      zeros(2,6) zeros(2,2)];
Q = [B1*B1'     zeros(6,2);
      zeros(2,6) zeros(2,2)]; 
R = blkdiag(C'*C,D'*D); 

%%%%optimization in the parameter space
W1 = blkdiag(sdpvar(3,3),sdpvar(3,3));
W2 = blkdiag(sdpvar(3,1),sdpvar(3,1));
W3 = sdpvar(2,2);
W = [W1  W2;
     W2' W3];
V=[eye(6)  zeros(6,2)];
z=trace(R*W);
MyLMI = [W >= 0, V*(F*W + W*F'+Q)*V' <= 0];
options=sdpsettings('solver','sdpt3');%sdpt3
optimize(MyLMI, z,options)
J=value(z)
K=value(W2)'*value(W1)^-1
% eig(A-B2*K)