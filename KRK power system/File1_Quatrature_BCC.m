clear all;  clear memory; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boyd/Clenshaw-Curtis Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  f=10; Lq=2.43;
% den=acot((f*2*pi)/(Lq));
% num=pi-acot((f*2*pi)/(Lq));
% Nq=(num) / (den)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ajuste1
w_selected=10;        % numero de puntos de 0 a pi Np=Nq
L1=18.00;             % constante de truncamiento del dominio para P     U<V
L2=18.44;             % constante de truncamiento del dominio para Q
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ajuste2
% w_selected=20;       % numero de puntos de 0 a pi Np=Nq
% L1=9.00;             % constante de truncamiento del dominio para P     U<V
% L2=9.47;             % constante de truncamiento del dominio para Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ajuste3
% w_selected=30;       % numero de puntos de 0 a pi Np=Nq
% L1=6.00;             % constante de truncamiento del dominio para P     U<V
% L2=6.38;             % constante de truncamiento del dominio para Q
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %ajuste4
% w_selected=40;       % numero de puntos de 0 a pi Np=Nq
% L1=4.00;             % constante de truncamiento del dominio para P     U<V
% L2=4.82;             % constante de truncamiento del dominio para Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %ajuste5
% w_selected=50;       % numero de puntos de 0 a pi Np=Nq
% L1=3.00;             % constante de truncamiento del dominio para P     U<V
% L2=3.87;             % constante de truncamiento del dominio para Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ajuste6
%% ES EL DE MEJORES RESULTADOS
% w_selected=60;       % numero de puntos de 0 a pi Np=Nq
% L1=3.00;             % constante de truncamiento del dominio para P     U<V
% L2=3.23;             % constante de truncamiento del dominio para Q
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Np=w_selected;
Nq=w_selected;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_1          =(pi/(w_selected+1));         %espacio entre puntos (h)
t_l_1        =1*h_1:h_1:(w_selected)*h_1;  %nodos en cuadratura BCC (tao_j)
omega_1      =L1.*cot(t_l_1);              %puntos frecuenciales desordenados (omega_j)
weight_w_l_1= sqrt(((L1*pi)./((w_selected+1).*(sin(t_l_1).^2)))); %pesos desordenados (phi_j)
%ORDERING
for j=1:Np
    if j<=(Np/2)
       omega_j  (Np-2*j+1)=omega_1(j);
       phi_j    (Np-2*j+1)=weight_w_l_1(j);
    end
    if j>(Np/2)
       omega_j   (2*j-Np)=omega_1(j);
       phi_j     (2*j-Np)=weight_w_l_1(j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v
h_2         =(pi/(w_selected+1));          %espacio entre puntos
t_l_2       =1*h_2:h_2:(w_selected)*h_2;   %nodos en cuadratura BCC (tao_k)
zeta_1      =L2.*cot(t_l_2);               %puntos frecuenciales desordenados (zeta_k)
weight_w_l_2=sqrt(((L2*pi)./((w_selected+1).*(sin(t_l_2).^2)))); %pesos desordenados (rho_k)
%ORDERING
for k=1:Nq
    if k<=(Nq/2)
       zeta_k  (Nq-2*k+1)=zeta_1(k);
       rho_k   (Nq-2*k+1)=weight_w_l_2(k);
    end
    if k>(Nq/2)
       zeta_k   (2*k-Nq)=zeta_1(k);
       rho_k    (2*k-Nq)=weight_w_l_2(k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 fmin=omega_j(1,1)/(2*pi)
 fmax=zeta_k(1,end-1)/(2*pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save DATA_KRK omega_j zeta_k phi_j rho_k   %% GUARDO EL AJUSTE PERO DESPUES AGREGO MANUAL EN EL Fichero 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%