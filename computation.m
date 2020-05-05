A=-2;
B=1;
C=[4;-5];
D=[0;0];
Bf=1;
Df=[6;1];
Bd=[4,-4];
Dd=[-4,0;3,2];
[U,Sigma,V]=svd(Df);
U1=U(:,1);
U2=U(:,2);
S1=Sigma(1,:);
Df1=V*U1'/S1;
Df2=U2';
A1=A-Bf*Df1*C;
C2=Df2*C;
B1=Bd-Bf*Df1*Dd;
D2=Df2*Dd;
C1=Df1*C;
C2=Df2*C;
D1=Df1*Dd;
D2=Df2*Dd;

%CVX
cvx_begin
variable P(1,1) symmetric
variable Z(1,1)
variable S(1,1)
variable g2
m4=[A1*P+P*A1+Z*C2+C2*Z,P*B1(1)+Z*D2(1),P*B1(2)+Z*D2(2),C1+C2*S;
    B1(1)*P+D2(1)*Z,-g2,0,D1(1)+D2(1)*S;
   B1(2)*P+D2(2)*Z,0,-g2,D1(2)+D2(2)*S;
   C1+S*C2,D1(1)+S*D2(1),D1(2)+S*D2(2),-1];
minimize (g2)
subject to
-m4==semidefinite(4);
cvx_end

R=Z/P;
H=Df1+S*Df2;
L=-Bf*Df1+R*Df2;
g=sqrt(g2);
if (A+L*C)>0
    X=2*(A+L*C)/((H*C)^2);
    L=L-X*(H*C)'*H;
end
A_F=A+L*C;
C_F=H*C;

% LMI solver
% setlmis([]);
% g2=lmivar(1,[1,1]);P=lmivar(1,[1,1]);
% Z=lmivar(2,[1,1]);S=lmivar(2,[1,1]);
% L=newlmi;
% lmiterm([L,1,1,P],1,A1,'s');lmiterm([L,1,1,Z],1,C2,'s');
% lmiterm([L,1,2,P],1,B1);lmiterm([L,1,2,Z],1,D2);
% lmiterm([L,1,3,0],C1');lmiterm([L,1,3,-S],C2',1);
% lmiterm([L,2,2,g2],-.5,1,'s');
% lmiterm([L,2,3,0],D1');lmiterm([L,2,3,-S],D2',1);
% lmiterm([L,3,3,0],-1);
% LMI=getlmis;
% [g2,x]=mincx(LMI,eye(decnbr(LMI),1));
% P=dec2mat(LMI,x,P);
% Z=dec2mat(LMI,x,Z);S=dec2mat(LMI,x,S);
% R=Z/P;
% H=Df1+S*Df2;
% L=-Bf*Df1+R*Df2;
% L=L-2*(A+L*C)*H/(H*C);
% A_F=A+L*C;
% C_F=H*C;