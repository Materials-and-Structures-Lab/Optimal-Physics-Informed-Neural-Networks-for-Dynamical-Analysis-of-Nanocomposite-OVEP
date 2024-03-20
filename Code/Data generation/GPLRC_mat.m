format long
global  A11 A12 A22 A66 B11 B12 B22 B66 D11 D12 D22 D66 rho  
global  a b h phiNx1 phiNy1 phiMx1 phiMy1 

syms z y x
%%%%% Temperature
T = 300;
deltaT = T - 300;
%%%%% Geometrical parameters
h = 0.005;
a = 0.3;
b = a; 
NL = 12;    

%%%%% The material properties of matrix and Graphene Platelets
Em = (4854.6 - 6.1816*deltaT)*10^9;
EG = (1087.8 - 0.261*deltaT)*10^9;
rhom = 1200;
rhoG = 1062.5;
alpham = 60*10^(-6);
alphaG = (13.92 - 0.0299*T) * 10^(-6);
vm = 0.34;
vG = 0.186;
%%%%%%%%%%%%%%%%%%%%%
aG = 2.5 * (10^(-6));
bG = 1.5 * (10^(-6));
tG = 1.5 * (10^(-9));
%%%%%%%%%%%%%%%%%%%%%
eL = 2*aG/tG;
eT = 2*bG/tG;
nL = (EG/Em-1)/(EG/Em+eL);
nT = (EG/Em-1)/(EG/Em+eT);
%%%%% Stiffness 
A11 =  0;
A12 =  0;
A22 =  0;
A66 =  0;
%%%%%%%%%%%%%%%%%%%%%
B11 =  0;
B12 =  0;
B22 =  0;
B66 =  0;
%%%%%%%%%%%%%%%%%%%%%
D11 =  0;
D12 =  0;
D22 =  0;
D66 =  0;
%%%%%%%%%%%%%%%%%%%%%
rho = 0;
%%%%%%%%%%%%%%%%%%%%%
phiNx1 =  0;
phiNy1 =  0;
phiMx1 =  0;
phiMy1 =  0;
%%%%% Layer Divide
e =  h/NL;
u = (-h/2) - e;
o = -h/2;
%%%%% ABD for laminated plate
for k = 1:1:NL
u = u + e;
o = o + e;
WG = 0.01; %%% GPL's Weight
%%%%%%%%%%%%%% GPL distribution   
%%%%%%%%%%%%%% UD
% WGk = WG;
%%%%%%%%%%%%%% FG-O
% WGk=4* WG *(0.5 + abs(k - 0.5*(NL + 1))) / (2 + NL);
%%%%%%%%%%%%%% FG-X
WGk=4* WG *((NL/2 + 0.5) - abs(k - (NL/2 + 0.5))) / (2 + NL);
%%%%%%%%%%%%%% Volume fraction of GPL
VGk = WGk / (WGk + (rhoG/rhom) * (1-WGk));
%%%%%%%%%%%%%% Material properties of kth layer
Ek = (3/8).*Em.*(1+(-1).*nL.*VGk).^(-1).*(1+eL.*nL.*VGk)+(5/8).*Em.*(1+(-1).* ...
  nT.*VGk).^(-1).*(1+eT.*nT.*VGk);
rhok = @(z) rhom.*(1+(-1).*VGk)+rhoG.*VGk;
vk = vG.*VGk+(1+(-1).*VGk).*vm;
alphak = alpham.*(1+(-1).*VGk)+alphaG.*VGk;
%%%%%%display ket qua
rho111 = rhom.*(1+(-1).*VGk)+rhoG.*VGk;
% disp(Ek)
% disp(rho111)
% disp(vk)
%%%%%%%%%%%%%% Elastic stiffness Components
Q11a = @(z) Ek.*(1+(-1).*vk.^2).^(-1);
Q12a = @(z) Ek.*vk.*(1+(-1).*vk.^2).^(-1);
Q22a = @(z) Ek.*(1+(-1).*vk.^2).^(-1);
Q66a = @(z) (1/2).*Ek.*(1+vk).^(-1);

Q11b = @(z) Ek.*(1+(-1).*vk.^2).^(-1).*z;
Q12b = @(z) Ek.*vk.*(1+(-1).*vk.^2).^(-1).*z;
Q22b = @(z) Ek.*(1+(-1).*vk.^2).^(-1).*z;
Q66b = @(z) (1/2).*Ek.*(1+vk).^(-1).*z;

Q11d = @(z) Ek.*(1+(-1).*vk.^2).^(-1).*z.^2;
Q12d = @(z) Ek.*vk.*(1+(-1).*vk.^2).^(-1).*z.^2;
Q22d = @(z) Ek.*(1+(-1).*vk.^2).^(-1).*z.^2;
Q66d = @(z) (1/2).*Ek.*(1+vk).^(-1).*z.^2;
%%%%%%%%%%%%%%% ABD Matrix
A11k = integral(Q11a,u,o,'ArrayValued',true);
A11 = A11 + A11k;
A12k = integral(Q12a,u,o,'ArrayValued',true);
A12 = A12 + A12k;
A22k = integral(Q22a,u,o,'ArrayValued',true);
A22 = A22 + A22k;
A66k = integral(Q66a,u,o,'ArrayValued',true);
A66 = A66 + A66k;

D11k = integral(Q11d,u,o,'ArrayValued',true);
D11 = D11 + D11k;
D12k = integral(Q12d,u,o,'ArrayValued',true);
D12 = D12 + D12k;
D22k = integral(Q22d,u,o,'ArrayValued',true);
D22 = D22 + D22k;
D66k = integral(Q66d,u,o,'ArrayValued',true);
D66 = D66 + D66k;
%%%%%%%%%%%%%%% Density
rho1a = integral(rhok,u,o,'ArrayValued',true);
rho = rho + rho1a;
%%%%%%%%%%%%%%% Temperature factor
% Q11temper = @(z) (Q11a(z)*alphak + Q12a(z)*alphak)*deltaT;
% Q12temper = @(z) (Q12a(z)*alphak + Q22a(z)*alphak)*deltaT;
Q11temper = @(z) (Q11a(z)*alphak + Q12a(z)*alphak)*deltaT;
Q12temper = @(z) (Q12a(z)*alphak + Q22a(z)*alphak)*deltaT;
phiNxk = integral(Q11temper,u,o,'ArrayValued',true);
phiNx1 = phiNx1 + phiNxk;
phiNyk = integral(Q12temper,u,o,'ArrayValued',true);
phiNy1 = phiNy1 + phiNyk;
end