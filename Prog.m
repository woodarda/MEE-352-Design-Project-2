close all; clear; clc;

%% Chem Properties %%

Ti= 80+273; % Inlet Temperature (K)
To= 40+273; % Outlet Temperature (K)
mdot1=3200/(90*60); % Chemical flow rate Scenario 1 (kg/s)
mdot2=3200/(150*60); % Chemical flow rate Scenario 2 (kg/s)

[Temp, P,vf, hfg, cpf, muf, kf, Prf]=AW_Interpolation((To+Ti)*0.5);

Chem=table(Ti,To,mdot1,mdot2,Temp,P,vf,hfg,cpf,muf,kf,Prf);

%% Water Properties %%
Ti=10+273; % water inlet temperature
mdot=0.5; % Water flow rate (kg/s)
To=Chem.To;
Water=table(Ti, To, mdot);

qh=[Chem.mdot1,Chem.mdot2].*Chem.cpf*(Chem.Ti-Chem.To);


Water.To=CalcToH2O(Water.To,Water.Ti,Water.mdot,Chem.cpf,Chem.mdot1);

Tbara=Water.To/2+Water.Ti/2;

[Temp, Water.P,Water.vf, Water.hfg, Water.cpf, Water.muf, Water.kf, ...
    Water.Prf]=AW_Interpolation(Tbara);

Chem.Ch=Chem.mdot1*Chem.cpf;
Water.Cc=Water.mdot*Water.cpf;

