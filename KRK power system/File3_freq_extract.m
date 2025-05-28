clear all; close all; clear memory; clc;
d16m_mpm;                              %% Cargo la data
File2_quadrature_sorting;        %% abro los puntos frecuenciales y los pesos
w_k_temporal=[omega_jota';zeta_ka'];   %% los acomodo para disminuir las simulaciones y hacerlo en una corrida
w_k=reshape(w_k_temporal,1,60);        %%
w_k=w_k/(2*pi); % Convirtiendo a Hertz (Vector de frecuencias)

global maqui_num frec ttt amplitud fase Mag_bus
[mm,nn]= size(mac_con);
Numdbus= length(bus(:,1));
amplitud=zeros(mm,Numdbus+mm,length(w_k));
fase=zeros(mm,Numdbus+mm,length(w_k));

cont1=1;cont5=1;
for frec = w_k
    for maqui_num=1:16  %aqui pongo manual el numero de generadores que tienen sistema de exitacion porque en otros sistemas mas complejos modelados en past como el NPCC no todos lo tienen
        [pm_sig mac_spd Mag_bus t]=s_simus();
        dt=t(2)-t(1);Fs=1/dt;
        sigau=pm_sig(cont1,122:end);
        sigau=sigau-mean(sigau);
        Nm=length(sigau);
        
        [Magnitud, Fase]=fft_ex(Fs,sigau,Nm);
%Magnitud
%Fase
        %%
        cont2=1;cont3=1;cont4=1;
        for Aux1 = 1:(Numdbus+mm)
            if Aux1 <= mm
                sigaux=mac_spd(cont2,122:end);
                sigaux=sigaux-mean(sigaux);
                Nm=length(sigaux);
                [Mag2, Fase2]=fft_ex(Fs,sigaux,Nm);
%Mag2
%Fase2
                AMPLI=Mag2/Magnitud;
                FASE=Fase2-Fase;
                amplitud(cont1,cont3,cont5)=AMPLI;
                fase(cont1,cont3,cont5)=FASE(1,1);
                cont2=cont2+1;
                cont3=cont3+1;
            else
                 sigaux=Mag_bus(cont4,122:end);
                sigaux=sigaux-mean(sigaux);
                Nm=length(sigaux);
                [Mag2, Fase2]=fft_ex(Fs,sigaux,Nm);
%Mag2
%Fase2
                AMPLI=Mag2/Magnitud;
                FASE=Fase2-Fase;
                amplitud(cont1,cont3,cont5)=AMPLI;
                fase(cont1,cont3,cont5)=FASE(1,1);
                cont4=cont4+1;
                cont3=cont3+1;
            end
        end
        cont1=cont1+1;
    end
    cont5=cont5+1;
    cont1=1;
end
save DATA_NETS_amp_fase_wk_hz amplitud fase w_k