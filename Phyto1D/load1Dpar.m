function par = load1Dpar()
%defines most parameters for Phyto_1D model
%themodynamic parameters
par.TC = 25;                 %Temp in degrees C
par.T = par.TC + 273.15;         %Temp in degrees K
par.S = 35;                 %Salinity in psu
par.R = 8.314E-3;           % gas constant in kJ/(mol K);
DCF = diff_coef(par.TC,par.S, 1);        %diffusion coefficients (in cm2/s) at current temp, salinity, final variable is pressure: assume 1 atm.
DCF = DCF ./ 100; %convert from cm2/s to dm2/sd
par.Dco2  = DCF(3);     %diffusivity of CO2 in dm2/s
par.Db  = DCF(8);       %diffusivity of HCO3- in dm2/s
par.Dc  = DCF(9);       %diffusivity of CO32- in dm2/s

%phyto parameters
par.Rc = 10/1E5;     %radius of phytoplankton; in um converted to dm;
par.Vc = (4/3) * pi * par.Rc^3;
par.DBL = 100/1E5; %length of diffusive boundary layer in um converted to dm; a few times the radius of the phytoplankton ;10x cell radius seems to by fine, numerical problems arise around 30x cell radius, and become major by 100x 
par.shells(1) = 30;     %number of shells to split cell into 
par.shells(2) = 75;     %number of shells to split DBL into
par.L = par.Rc + par.DBL;
par.st(1) = par.Rc ./ par.shells(1);    %shell thickness in cell
par.st(2) = par.DBL ./par.shells(2);    %shell thickness in DBL
par.stInt = (par.st(1) + par.st(2))./2; %interface distance - between center of innermost shell in DBL and outermost shell in cell (the membrane)

%equilibrium constants for DIC, water, and borate
par.K1 = 10^-((3633.86./par.T) - 61.2172 + 9.6777*log(par.T) - 0.011555 .* par.S + 0.0001152 .* (par.S^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
par.K2 = 10^((-471.78/par.T) - 25.929 + 3.16967 * log(par.T) + 0.01781 * par.S - 0.0001122*par.S^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
par.Kw = exp(148.96502 - (13847.26 ./ par.T) - 23.6521 .* log(par.T) + (par.S^.5).*((118.67 ./ par.T) - 5.977 + 1.0495 .* log(par.T)) - 0.01615 .* par.S); %ion product of water, CO2 methods DOE 1994
par.Kb = exp(((-8966.9 - 2890.53 .* par.S^0.5 - 77.942 .* par.S + 1.728 .* par.S^1.5 - 0.0996 .* par.S^2)./par.T) + 148.0248 + 137.1942 .* par.S^0.5 + 1.62142 .* par.S...
    -(24.4344 + 25.085 .* par.S^0.5 + 0.2474 .* par.S) .* log(par.T) + 0.053105 .* par.S^0.5 .* par.T); %Kb boric acid/borate equilibrium constant

%kinetic parameters
par.kp1 = exp(1246.98 - 6.19E4 ./ par.T - 183 .* log(par.T));         % CO2 + H2O -> H+ + HCO3- Johnson 1982; as presented in Zeebe and Wolf-Gladrow
par.kp4 = 4.70E7 .* exp(-23.2 ./ (8.314E-3 .* par.T));                % CO2 + OH- -> HCO3- Johnson 1982
par.km1 = par.kp1./par.K1;                                            % H+ + HCO3- -> CO2 + H2O
par.km4 = par.kp4.*par.Kw./par.K1;                                    % HCO3- -> CO2  + OH-
par.kph5 = 5.0E10;                                                    % CO32- + H+ -> HCO3-; /(M s); Zeebe and Wolf Gladrow
par.kmh5 = par.kph5 .* par.K2;                                        % HCO3- -> CO32- + H+;
par.kpoh5 = 6.0E9;                                                    % HCO3- + OH- -> CO32- + H2O; Zeebe and Wolf-Gladrow
par.kmoh5 = par.kpoh5 .* par.Kw ./ par.K2;                            % CO32- + H2O -> HCO3- + OH-;

%solution parameters, pH which is constant, and initial DIC system
%conditions
par.H   = 10^-8.046;        %H+ concentration
par.OH  = (par.Kw./par.H);        %OH- concentration
par.CO2 = 11.243E-6;         %CO2 molar  11.24E-6
par.B   = 1775.5E-6;        %HCO3- molar
par.C   = 213.3E-6;         %CO32- molar 213.3E-6