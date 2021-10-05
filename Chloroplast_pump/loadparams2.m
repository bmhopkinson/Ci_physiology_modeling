function p = loadparams2(p,pFile)


fp = fopen(pFile,'r');

if fp==-1
   error(['File ' pFile ' not found or permission denied.']);
end

i=1;


while 1 
    line = fgetl(fp);
    if (strcmp('________',line))
        break
    end
    raw{i} = line;
    i = i+1;
end

for i = 1:length(raw)
    [id] = sscanf(raw{i},'%s[\t]%*s');
    
    switch id
        
        case 'pHc'
            p.pHc = sscanf(raw{i},'%*s %f');        %cytoplasmic pH
        case 'pHp'
            p.pHp = sscanf(raw{i},'%*s %f');        %chloroplast stroma pH
        case 'kcf'
            p.kcf = sscanf(raw{i},'%*s %e');        %cytoplasmic CA activity
        case 'kpf'    
            p.kpf = sscanf(raw{i},'%*s %e');        %chloroplast stroma CA activity
        case 'kyf'    
            p.kyf = sscanf(raw{i},'%*s %e');        %pyrenoid CA activity
        case 'mRub'    
            p.mRub = sscanf(raw{i},'%*s %e');        %mols of RubisCO molecules
        case 'Vm_Bc'    
            p.Vm_Bc = sscanf(raw{i},'%*s %e');        %maximum uptake rate of cytoplasmic HCO3- transporter in mol/cell/s
        case 'Km_Bc'    
            p.Km_Bc = sscanf(raw{i},'%*s %e');        %Km for cytoplasmic HCO3- transporter (mol/cm3)
        case 'Vm_Bp'    
            p.Vm_Bp = sscanf(raw{i},'%*s %e');        %maximum uptake rate of chloroplast pump HCO3- transporter in mol/cell/s
        case 'Km_Bp'    
            p.Km_Bp = sscanf(raw{i},'%*s %e');        %Km for chloroplast pump HCO3- transporter (mol/cm3)
        case 'fc_c'    
            p.fc_c = sscanf(raw{i},'%*s %e');        %MTC for CO2 into cytoplasm, from Hopkinson et al. 2011
        case 'fb_c'    
            p.fb_c = sscanf(raw{i},'%*s %e');        %MTC for HCO3- into cytoplasm, generally set to zero
        case 'fc_p'    
            p.fc_p = sscanf(raw{i},'%*s %e');        %MTC for CO2 into chloroplast, from Hopkinson et al. 2011
        case 'fb_p'    
            p.fb_p = sscanf(raw{i},'%*s %e');        %MTC for HCO3- into chloroplast, generally set to zero
        case 'fc_y'    
            p.fc_y = sscanf(raw{i},'%*s %e');        %MTC for CO2 into pyrenoid, from Hopkinson et al. 2011
    end
end


%define pameters for CCM mechanistic model
p.T  = p.TC + 273.15;
p.S  = 35;          %salinity
p.H  = 10^-p.pH;    %H+ of extracellular solution  
p.Hc  = 10^-p.pHc;  %H+ concentration in the cytoplasm;
p.Hp = 10^-p.pHp;   %H+ concentration in the stroma/pyrenoid

%equilibrium constants for DIC, water, and borate
p.K1 = 10^-((3633.86./p.T) - 61.2172 + 9.6777*log(p.T) - 0.011555 .* p.S + 0.0001152 .* (p.S^2));           %K1 equilibrium constant between CO2 and HCO3; Lueker, Dickson, Keeling Mar Chem. 2000.
p.K2 = 10^((-471.78/p.T) - 25.929 + 3.16967 * log(p.T) + 0.01781 * p.S - 0.0001122*p.S^2);         %K2 equilbrium constant from Lueker, Dickson, Keeling Mar Chem 2000
p.Kw = exp(148.96502 - (13847.26 ./ p.T) - 23.6521 .* log(p.T) + (p.S^.5).*((118.67 ./ p.T) - 5.977 + 1.0495 .* log(p.T)) - 0.01615 .* p.S); %ion product of water, CO2 methods DOE 1994
p.Kb = exp(((-8966.9 - 2890.53 .* p.S^0.5 - 77.942 .* p.S + 1.728 .* p.S^1.5 - 0.0996 .* p.S^2)./p.T) + 148.0248 + 137.1942 .* p.S^0.5 + 1.62142 .* p.S...
    -(24.4344 + 25.085 .* p.S^0.5 + 0.2474 .* p.S) .* log(p.T) + 0.053105 .* p.S^0.5 .* p.T); %Kb boric acid/borate equilibrium constant
p.bfrac_e = (1/ (1+ p.K2./p.H));       % fraction of "B" pool that is HCO3- in extracellular solution
p.bfrac_i = (1/ (1+ p.K2./p.Hc));       % fraction of "B" pool that is HCO3- in cytoplasm
p.bfrac_x = (1/ (1+ p.K2./p.Hp));       % fraction of "B" pool that is HCO3-   in chloroplast and pyrenoid

%kinetic pameters
p.kp1 = exp(1246.98 - 6.19E4 ./ p.T - 183 .* log(p.T));         % CO2 + H2O -> H+ + HCO3- Johnson 1982; as presented in Zeebe and Wolf-Gladrow
p.kp4 = 4.70E7 .* exp(-23.2 ./ (8.314E-3 .* p.T));                % CO2 + OH- -> HCO3- Johnson 1982
p.km1 = p.kp1./p.K1;                                            % H+ + HCO3- -> CO2 + H2O
p.km4 = p.kp4.*p.Kw./p.K1;                                    % HCO3- -> CO2  + OH-
p.kph5 = 5.0E10;                                                    % CO32- + H+ -> HCO3-; /(M s); Zeebe and Wolf Gladrow
p.kmh5 = p.kph5 .* p.K2;                                        % HCO3- -> CO32- + H+;
p.kpoh5 = 6.0E9;                                                    % HCO3- + OH- -> CO32- + H2O; Zeebe and Wolf-Gladrow
p.kmoh5 = p.kpoh5 .* p.Kw ./ p.K2;                            % CO32- + H2O -> HCO3- + OH-;

%diffusion coefficients
DCF = diff_coef(p.TC,p.S, 1);        %diffusion coefficients (in cm2/s) at current temp, salinity, final variable is pressure: assume 1 atm.
p.Db  = DCF(8);       %diffusivity of HCO3- in cm2/s

%CO2 hydration/HCO3- dehydration rates 
p.kuf = p.kp1 + p.kp4.*(p.Kw./p.H);     %CO2 hydration rate in bulk solution
p.kur = p.bfrac_e.* p.kuf .* (p.H./p.K1);           %HCO3- dehyration rate in bulk solution
p.kcr = p.bfrac_i.* p.kcf .* (p.Hc./p.K1);           %HCO3- dehyration rate in cytoplasm (/s)
p.kpr = p.bfrac_x.* p.kpf .* (p.Hp./p.K1);           %HCO3- dehyration rate in chloroplast stroma (/s)
p.kyr = p.bfrac_x.* p.kyf .* (p.Hp./p.K1);           %HCO3- dehyration rate in pyrenoid (/s)

%volumes
p.Ve = 1;                               %volume of solution (cm3);
p.Vc = 6.6E-11;                         %total volume of a single cell (cm3);
p.Vp = 0.6E-11;                         %chloroplast stroma volume
Rpyr = 0.3/1E4;                         %radius of the pyrenoid based on BCA localization, Hopkinson et al. 2011
Npyr = 1.9;                             %number of pyrenoids per cell, Hopkinson et al. 2011
p.Vy = Npyr.*(4.*pi./3).*(Rpyr^3);        %volume of the pyrenoid (cm3), as estimated in Hopkinson et al 2011 from BCA clusters;

%enzyme kinetics
p.kcat_R = 3.4;                   %Pt Rubisco turnover rate (/s) according to Whitney et al. 2001 Plant Journal 26: 535-547
p.Km_R   = 41 .* 1E-9;            %Pt Rubisco Km (uM converted to mol/cm3) according to Badger et al. 1998

%additional mass transfer coefficients
p.fb_y = Npyr.*4.*pi.*p.Db.*Rpyr;                     %MTC for HCO3- into pyrenoid.

end


