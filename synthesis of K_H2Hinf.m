m1=0.128;
m2=1.977;
b1=17.7;
b2=19.6;
k1=6797;
k11=k1*1.2;
k12=k1*0.8;
kf1=17.7;
kv1=0.2;
kf2=24;
kv2=1;
f=4.8;
A1=[0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 -k11/m1 -b1/m1 0 k11/m1 b1/m1;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   0 k11/m2 b1/m2 0 -k11/m2 -b2/m2-b1/m2];%之前-b2/m2-b1/m2推错了
A2=[0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 -k12/m1 -b1/m1 0 k12/m1 b1/m1;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   0 k12/m2 b1/m2 0 -k12/m2 -b2/m2-b1/m2];
B2=[0 0;
    0 0;
    -kf1*kv1/m1 0;
    0 0;
    0 0;
    kf1*kv1/m2 -kf2*kv2/m2];
B1=[0;0;0;0;0;f/m2];%eye(6);
%   B1=eye(6);
C=[100000 0 0 0 0 0;
   0 0 0 1000 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0];
D=[0 0 ;
   0 0 ;
   1 0 ;
   0 1 ;];
F1 = [A1          -B2;
      zeros(2,6) zeros(2,2)];
  F2 = [A2          -B2;
      zeros(2,6) zeros(2,2)];
Q = [B1*B1'     zeros(6,2);
      zeros(2,6) zeros(2,2)]; 
R = blkdiag(C'*C,D'*D); 
W1 = blkdiag(sdpvar(3,3),sdpvar(3,3));
W2 = blkdiag(sdpvar(1,3),sdpvar(1,3))';
W3 = sdpvar(2,2);
W = [W1  W2;
     W2' W3];
 V=[eye(6)  zeros(6,2)];
gamma=0.25;
z=gamma^2*trace(R*W);
% MyLMI = [W >= 0,  V*(F*W+W*F'+W*R*W+gamma^-2*Q)*V'<= 0];
MyLMI = [W >= 0, [eye(8)        R^(1/2)*W*V';
                  V*W*R^(1/2)    -V*F1*W*V'-V*W*F1'*V'-gamma^-2*V*Q*V'] >= 0,...
                  [eye(8)        R^(1/2)*W*V';
                  V*W*R^(1/2)    -V*F2*W*V'-V*W*F2'*V'-gamma^-2*V*Q*V'] >= 0 
                  ];
% MyLMI = [W >= 0, [-V*F*W*V'-V*W*F'*V'-gamma^-2*V*Q*V'  V*W*R^(1/2);
%                    R^(1/2)*W*V'   eye(8)] >= 0];
options=sdpsettings('solver','sdpt3');
optimize(MyLMI, z,options)
J=value(z)
K=value(W2)'*value(W1)^-1
%eig(A-B2*K)
sqrt(J)
norm(ss(A1-B2*K,B1,C-D*K,zeros(4,1)),2)
norm(ss(A2-B2*K,B1,C-D*K,zeros(4,1)),2)
% [norm(ss(A1-B21*K,B1,C-D*K,zeros(4,1)),2),...
%  norm(ss(A1-B22*K,B1,C-D*K,zeros(4,1)),2),...
%  norm(ss(A1-B23*K,B1,C-D*K,zeros(4,1)),2),...
%  norm(ss(A1-B24*K,B1,C-D*K,zeros(4,1)),2),...
%  norm(ss(A2-B21*K,B1,C-D*K,zeros(4,1)),2),...
%  norm(ss(A2-B22*K,B1,C-D*K,zeros(4,1)),2),...
%  norm(ss(A2-B23*K,B1,C-D*K,zeros(4,1)),2)...
%  norm(ss(A2-B24*K,B1,C-D*K,zeros(4,1)),2)],
% K=[-99.9997615151877,-0.273888173701886,-0.260075534308896,0,0,0;0,0,0,-100.021904980375,-799.529956282572,-10.9491391697054];
%  K=[-99999.9927572828,-374.876540055279,-2.90050224538452,0,0,0;0,0,0,-1000.09362251067,-2154.34846613427,-20.1339782524757];
[norm(ss(A1-B2*K,B1,C-D*K,zeros(4,1)),inf),...
 norm(ss(A2-B2*K,B1,C-D*K,zeros(4,1)),inf)]