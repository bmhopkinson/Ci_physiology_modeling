function ydot = D1deriv(time, y, par)
% derivative function for kinetic description of the CO2 system in
% seawater, assumes pH is constant. the only important kinetic terms are
% hydration and dehydration of CO2/HCO3-


tshells = (par.shells(1) + par.shells(2));
P = 1E-15; %photosynthetic rate in mol/cell/s
Pv = P ./ ((4*pi/3) * (par.Rc - par.st(1))^3); %volumentric photosynthetic rate mol/L/s NOTE that photosyntheis is negeletcted in the membrane layer so volume used is only the intracellular space
Pv = [Pv; 0; 0]; %put in vector form

%formulate reaction matrix in DBL (outside cell)
Ro(1,1:3) = [-par.kp1 - par.kp4 .* par.OH, par.km1 .* par.H + par.km4, 0];       %for dCO2/dt
Ro(2,1:3) = [ par.kp1 + par.kp4 .* par.OH, -par.km1 .* par.H - par.km4 - par.kmh5 - par.kpoh5 .* par.OH, ...
        par.kph5 .* par.H + par.kmoh5];                         % for dHCO3-/dt
Ro(3,1:3) = [ 0, par.kmh5 + par.kpoh5 .* par.OH, - par.kph5 .* par.H - par.kmoh5];   % for dCO32-/dt

%formulate reaction matrix inside cell
kcf = par.kp1 .* 1;      %simulate CA activity inside the cell
kcr = par.km1 .* 1;
Ri(1,1:3) = [-kcf - par.kp4 .* par.OH, kcr .* par.H + par.km4, 0];       %for dCO2/dt
Ri(2,1:3) = [ kcf + par.kp4 .* par.OH, -kcr .* par.H - par.km4 - par.kmh5 - par.kpoh5 .* par.OH, ...
        par.kph5 .* par.H + par.kmoh5];                         % for dHCO3-/dt
Ri(3,1:3) = [ 0, par.kmh5 + par.kpoh5 .* par.OH, - par.kph5 .* par.H - par.kmoh5];   % for dCO32-/dt

R0 = [0, 0, 0; 0, 0, 0; 0, 0, 0];       %stop all chemical reactions

%formulate diffusion matrix
dv = (4*pi) .* [par.Dco2; par.Db; par.Dc];  %diffusivity vector, include 4*pi term
D = diag(dv);

%diffusion matrix for membrane layer
Dmem = diag(4*pi .* [par.Dco2; 0 ; 0]);

%compute derivative vector in cell
for i = 1:par.shells(1)  
    
    R = Ri;             %chose correct reaction matrix
    j = 3*(i-1)+1;      %starting index
    k = 3*(i-1)+3;      %ending index, 3 species are treated (CO2, HCO3-, CO32-)

    %inner radius and outer radius of shell
    ri = par.st(1) .* (i-1);
    ro = par.st(1) .* i;
    
    %volume of shell (or sphere for inner-most layer)
    Vsh = (4*pi/3).*(ro^3 - ri^3);
 
     % compute derivative vector (d[C]/dt, M/s)
    if (i == 1)
         ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).* (D *    (ro^2 .* (y(j+3:k+3,1) - y(j:k,1))./par.st(1))) - Pv;  %inner-most sphere of cell
    elseif (i == par.shells(1)-1)  
        ydot(j:k,1)  = (R * y(j:k,1)) + (1/Vsh).* (Dmem * (ro^2 .* (y(j+3:k+3,1) - y(j:k,1))./par.st(1)) + D * (ri^2 .* (y(j-3:k-3,1) - y(j:k,1))./par.st(1))) - Pv;   %shell just within membrane - account of lack of HCO3-/CO32- diffusion into membrane
    elseif (i == par.shells(1))
         ydot(j:k,1) =                  (1/Vsh).* (Dmem * (ro^2 .* (y(j+3:k+3,1) - y(j:k,1))./par.stInt       + ri^2 .* (y(j-3:k-3,1) - y(j:k,1))./par.st(1)));      %membrane layer, %treat as the cell membrane, permeable to CO2, but not to HCO3- or CO32-, neglect rxns in the membrane b/c rates and Ci concentrations (other than CO2)are ill-defined
    else
        ydot(j:k,1)  = (R * y(j:k,1)) + (1/Vsh).* (D *    (ro^2 .* (y(j+3:k+3,1) - y(j:k,1))./par.st(1)       + ri^2 .* (y(j-3:k-3,1) - y(j:k,1))./par.st(1))) - Pv;      %intermediate shells
    end
end

%compute derivative vector in DBL

for i = par.shells(1)+1:tshells
%for i = 1:par.shells(2)   
    R = Ro;             %chose correct reaction matrix
    j = 3*(i-1)+1;      %starting index
    k = 3*(i-1)+3;      %ending index, 3 species are treated (CO2, HCO3-, CO32-)
   
    %inner radius and outer radius of shell
    ri = par.Rc + par.st(2) .* (i - par.shells(1) - 1);
    ro = par.Rc + par.st(2) .* (i - par.shells(1));
    
    %volume of shell
    Vsh = (4*pi/3).*(ro^3 - ri^3);
    
    %eCA activity; 
    ksf = 0;     %forward rate constant in cm3/s
    ksr = ksf .* (par.H./par.K1);
    keCA = [-ksf, ksr, 0; ksf, -ksr, 0; 0, 0, 0]./1000;   %put rate constants in matrix form; and convert to dm3/s

     % compute derivative vector (d[C]/dt, M/s)
    if (i == par.shells(1)+1)
        ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).*(D * (ro^2 .* (y(j+3:k+3,1)     - y(j:k,1))./par.st(2)) + Dmem * (ri^2 .* (y(j-3:k-3,1) - y(j:k,1))./par.stInt)) + (keCA * y(j:k,1))./Vsh;  % shell against cell surface
    elseif (i == tshells)
        ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).*(D * (ro^2 .* (par.Cinit(1:3,1) - y(j:k,1))./par.st(2)          + ri^2 .* (y(j-3:k-3,1) - y(j:k,1))./par.st(2)));      %shell against bulk fluid
    else
        ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).*(D * (ro^2 .* (y(j+3:k+3,1)     - y(j:k,1))./par.st(2)          + ri^2 .* (y(j-3:k-3,1) - y(j:k,1))./par.st(2)));      %intermediate shells
    end
end
    
ydot;
    
return


