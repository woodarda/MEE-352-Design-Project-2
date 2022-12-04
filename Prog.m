close all; clear; clc;

%% To optimize Concentric tube alter these
Do=50/1000; % Inner tube's outer diameter
Di=45/1000; % Inner tube's Inner diameter
D=75/1000;  % Outer Tube's Diameter
ks=1000; % conduction coefficient W/m*K

%% to optimize Cross Flow Alter These
CcMixed=false;
Aligned=True;


% Dont Change Anything under here
%% Chem Properties %%

Ti= 80+273; % Inlet Temperature (K)
To= 40+273; % Outlet Temperature (K)
mdot1=3200/(90*60); % Chemical flow rate Scenario 1 (kg/s)
mdot2=3200/(150*60); % Chemical flow rate Scenario 2 (kg/s)



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

Run90CF=Run90;
Run150CF=Run150;

Run90CT=Run90;
Run150CT=Run150;
%% HEAT EXCHANGER: CONCENTRIC CIRCLES

[Run90CT,Run150CT]=ConcentricCircles(Run90CT,Run150CT,Chem,Water,Di,Do,D,ks);


%% HEAT EXCHANGER: CROSS FLOW
Run90CF.CcMixed=CcMixed;
Run150CF.CcMixed=CcMixed;

if Aligned==false
    Run90CF.Aligned=false;
    Run150CF.Aligned=false;
end

if Aligned==true
    Run90CF.Aligned=True;
    Run150CF.Aligned=True;
end

[Run90CF,Run150CF]=CrossFlow(Run90CF,Run150CF,Chem,Water,Di,Do,D,ks);
Run90CF
Run150CF

Run90CT
Run150CT