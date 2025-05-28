clear all; clear memory; close all; clc;
load('DATA_NETS_amp_fase_wk_hz.mat'); %fase amplitud w_k total ordenado y en Hz 
File2_cuadratura_ajuste_manual;       %Weights
phi_j=phi_jota';
rho_k=rho_ka';

[mm,nn,oo]=size(amplitud);
Hi=zeros(nn,mm,oo);
for cont1=1:oo
  for cont2=1:nn
      H=amplitud(:,cont2,cont1).*exp(fase(:,cont2,cont1).*1i);
      Hi(cont2,:,cont1)= H.';
  end
end
tic
si=0+(2*pi.*w_k.*1i);
si=si.';
Hs=Hi(:,:,:);  
s=si;
[m,p,N]=size(Hs);
cont1=1:2:N;
omega=s(cont1); 
cont2=2:2:N;
zeta =s(cont2);

for n=1:length(omega)
    H_P(:,:,n )=Hs(:,:,cont1(n));
end
for mm=1:length(zeta)
    H_Q(:,:,mm)=Hs(:,:,cont2(mm));
end

%% Aproximacion de P - H_P
omega1(1:2:2*length(omega),1)=omega;
omega1(2:2:2*length(omega),1)=conj(omega);
H_P_P(:,:,1:2:2*length(omega))=H_P;
H_P_P(:,:,2:2:2*length(omega))=conj(H_P);
omega=omega1;
H_P=H_P_P;


%% Aproximacion de Q - H_Q
zeta1(1:2:2*length(zeta),1)=zeta;
zeta1(2:2:2*length(zeta),1)=conj(zeta);
H_Q_Q(:,:,1:2:2*length(zeta))=H_Q;
H_Q_Q(:,:,2:2:2*length(zeta))=conj(H_Q);
zeta=zeta1;
H_Q=H_Q_Q;

omega=omega';
zeta=zeta';
L_complex=zeros(length(omega)*m,length(zeta)*p);
M_complex=zeros(length(omega)*m,length(zeta)*p);

outputs=m;
inputs=p;

for x=1:length(omega)
    for y=1:length(zeta)
           L_complex(1 + outputs*(x-1):outputs+outputs*(x-1),1 + inputs*(y-1):inputs+inputs*(y-1))=((H_P(:,:,x)-H_Q(:,:,y))/(omega(1,x)-zeta(1,y))).*(-1).*phi_j(1,x).*rho_k(1,y);
           M_complex(1 + outputs*(x-1):outputs+outputs*(x-1),1 + inputs*(y-1):inputs+inputs*(y-1))=((omega(1,x)*H_P(:,:,x)-zeta(1,y)*H_Q(:,:,y))/(omega(1,x)-zeta(1,y))).*(-1)*phi_j(1,x).*rho_k(1,y);   
           H_complex(outputs*(x-1) +1:outputs+outputs*(x-1),:)                                    =H_P(:,:,x).*phi_j(1,x);
           G_complex(:,inputs*(y-1)+1:inputs+inputs*(y-1))                                        =H_Q(:,:,y).*rho_k(1,y); 
    end
end
%% (L_complex, M_complex, H_complex, G_complex)
TL=double.empty;
TR=double.empty;
Ip_=eye(inputs); 
Ip=(1/sqrt(2)).*[Ip_ Ip_.*1j;Ip_ -Ip_.*1j];
Im_=eye(outputs); 
Im=(1/sqrt(2)).*[Im_ Im_;-Im_.*1j Im_.*1j];
for i = 1:length(omega)/2
    TL = blkdiag(TL,Im);  
end
for i = 1:length(zeta)/2
    TR = blkdiag(TR,Ip); 
end
L_real  = real(TL*L_complex*TR');             
M_real  = real(TL*M_complex*TR');
% h_real(:,:)  =L_real(:,end-inputs+1:end);
% gT_real(:,:) =M_real(end-outputs+1:end,:);
h_real  = real(TL*H_complex);
gT_real = real(G_complex*TR'); 
[Z,S,Y1] = svd(L_real,'econ');
Sv_ener = diag(S);
sv_ener=sum(Sv_ener);
En2=0.1;
flag2=0;
Ene2=0;
while En2<0.98     %% 98 con el ajuste de cuadratura utilizado arroja el mejor resultado para este sistema
    flag2=flag2+1;
    Ene2=Ene2+Sv_ener(flag2);
    En2=Ene2/sv_ener;
    ENE(flag2)=En2;
end
k=flag2;                    %% k singular values
%% REDUCTION
S1 = S(1:k,1:k);
Z1 = Z(:,1:k);
Z1=  Z1'; 
Y1= Y1(:,1:k);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S1=S1^(-0.5);
Er=     S1*Z1*L_real*Y1*S1;
Ar=     S1*Z1*M_real*Y1*S1;
Br=     S1*Z1*h_real;
Cr=gT_real*Y1 *S1;
toc
%% DEFINITION %%%%%%%%%%%%%%%%%%
Ar=inv(Er)*Ar;
Br=inv(Er)*Br;
Dr=zeros(size(Cr,1),size(Br,2));
Er=eye(size(Ar,1));
%%
sv_QBBT_29=diag(S)';

%%%% Frecuencuencies and damping ratios
TT=eig(Ar)
FREQUENCYS=imag(TT)/(2*pi)    %% frecuencies
DR=(-real(TT)./(abs(TT)));
DR=DR*100                     %% damping ratio


save QBBT_Model Ar Br Cr Dr Er
