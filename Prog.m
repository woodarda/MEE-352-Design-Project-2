close all; clear; clc;

%% Solid Properties %%

ks=1000; % conduction coefficient W/m*K


%% Chem Properties %%

Ti= 80+273; % Inlet Temperature (K)
To= 40+273; % Outlet Temperature (K)
mdot1=3200/(90*60); % Chemical flow rate Scenario 1 (kg/s)
mdot2=3200/(150*60); % Chemical flow rate Scenario 2 (kg/s)

Do=50/1000; % Inner tube's outer diameter
Di=25/1000; % Inner tube's Inner diameter
D=75/1000;  % Outer Tube's Diameter

% use interpolation to find the properties of the Temperature at Tbar
[Temp, P,v, hfg, cp, mu, k, Pr]=AW_Interpolation((To+Ti)*0.5);

% create a structure for the properties of the chemical
Chem=table(Ti,To,mdot1,mdot2,Temp,P,v,hfg,cp,mu,k,Pr);

%% Water Properties %%
Ti=10+273; % water inlet temperature
mdot=0.5; % Water flow rate (kg/s)
To=Chem.To; % initial guess for outlet temp of water
% create a structure for the properties of the water
Water=table(Ti, To, mdot);

% calculate the heat transfer to the interface between the two using
% chemical temp and properties
qh=[Chem.mdot1,Chem.mdot2].*Chem.cp*(Chem.Ti-Chem.To);

% iterate the water outlet temperature using CalcToH20
Water.To=CalcToH2O(Water.To,Water.Ti,Water.mdot,Chem.cp,Chem.mdot1);

% Find the actual average Temp of the water
Tbara=Water.To/2+Water.Ti/2;

% calculate the properties at the actual average temperature
[Temp, Water.P,Water.v, Water.hfg, Water.cp, Water.mu, Water.k, ...
    Water.Pr]=AW_Interpolation(Tbara);

% Caulculate Ch and Cc
Chem.Ch=[Chem.mdot1,Chem.mdot2].*Chem.cp;
Water.Cc=Water.mdot*Water.cp;

% Table for 90 min run
Run90=table(Chem.Ch(1),Water.Cc);
Run90=renamevars(Run90,["Var1","Var2"],["Ch","Cc"]);

% table for 2 75 min runs
Run150=table(Chem.Ch(2),Water.Cc);
Run150=renamevars(Run150,["Var1","Var2"],["Ch","Cc"]);

% Find Cmin and Cmax for 90 and 150 minute run
[Run90.Cmin,Run90.Cmax]=FindCminCmax(Run90.Ch,Run90.Cc);
[Run150.Cmin,Run150.Cmax]=FindCminCmax(Run150.Ch,Run150.Cc);

% heat exchange transfer for 90 and 150 run of the hot (Chemical) fluid
Run90.qh=Chem.cp*Chem.mdot1*(Chem.Ti-Chem.To);
Run150.qh=Chem.cp*Chem.mdot2*(Chem.Ti-Chem.To);

% Cr for the 90 min and 150 min run
Run90.Cr=Run90.Cmin/Run90.Cmax;
Run150.Cr=Run150.Cmin/Run150.Cmax;

% inlet Temp of the Chemical is the inlet hot temperature
Run90.Thi=Chem.Ti;
Run150.Thi=Chem.Ti;

% inlet Temp of the water is the inlet cold temperature
Run90.Tci=Water.Ti;
Run150.Tci=Water.Ti;

% Outlet temp of Chem is the hot outlet temperature 
Run90.Tho=Chem.To;
Run150.Tho=Chem.To;

% Outlet temp of Chem is the cold outlet temperature 
Run90.Tco=Water.To;
Run150.Tco=Water.To;

% Calculates the maximum potential heat exchange rate
Calc_qmax=@(Cmin,Thi,Tci) Cmin*(Thi-Tci);

% Run Calc_qmax for 90 and 150 minute scenarios
Run90.qmax=Calc_qmax(Run90.Cmin,Run90.Thi,Run90.Tci);
Run150.qmax=Calc_qmax(Run150.Cmin,Run150.Thi,Run150.Tci);

% Calculates the efficiency of the heat exhanger rate
CalcEff=@(qh,qmax) qh/qmax;

% Run CalcEff for 90 and 150 minute scenarios
Run90.Eff=CalcEff(Run90.qh,Run90.qmax);
Run150.Eff=CalcEff(Run150.qh,Run150.qmax);

% Calculates NTU if the flow is parallel
ParFlowNTU=@(Eff,Cr) -log(1-Eff*(1+Cr))/(1+Cr);


Run90.ParNTU=ParFlowNTU(Run90.Eff,Run90.Cr);
Run150.ParNTU=ParFlowNTU(Run150.Eff,Run150.Cr);

% Calculates NTU if the flow is Counter
CountFlowNTU=@(Eff,Cr) log((Eff-1)/(Eff*Cr-1))/(Cr-1);

% Runs CountNTU for both 90 and 150 scenarios
Run90.CountNTU=CountFlowNTU(Run90.Eff,Run90.Cr);
Run150.CountNTU=CountFlowNTU(Run150.Eff,Run150.Cr);

% Calculate Reynolds number, Nusselt number and convection coefficient 
% of the hot fluid
[Run90.Re_h,Run90.NuD_h,Run90.h_h]=IntFlowCCS(Chem,Chem.mdot1,Di);
[Run150.Re_h,Run150.NuD_h,Run150.h_h]=IntFlowCCS(Chem,Chem.mdot2,Di);

% Calculate Reynolds number, Nusselt number and convection coefficient 
% of the cold fluid
[Run90.Re_c,Run90.NuD_c,Run90.h_c]=IntFlowCCS(Water,Water.mdot,D-Do);
[Run150.Re_c,Run150.NuD_c,Run150.h_c]=IntFlowCCS(Water,Water.mdot,D-Do);

RConvFun=@(h,D) (h*pi*D)^-1;

Run90.Rconvh=RConvFun(Run90.h_h,Di);
Run90.Rcond=log(Do/Di)/(2*pi*ks);
Run90.Rconvc=RConvFun(Run90.h_c,Do);

R=[Run90.Rconvh,Run90.Rcond,Run90.Rconvc];

Run90.Req=sum(R);

Run150.Rconvh=RConvFun(Run150.h_h,Di);
Run150.Rcond=log(Do/Di)/(2*pi*ks);
Run150.Rconvc=RConvFun(Run150.h_c,Do);

R=[Run150.Rconvh,Run150.Rcond,Run150.Rconvc];

Run150.Req=sum(R);

