function ydot = D1Isoxderiv(time, y, par)
% derivative function for kinetic description of the CO2 system in
% seawater, assumes pH is constant. the only important kinetic terms are
% hydration and dehydration of CO2/HCO3-

tshells = (par.shells(1) + par.shells(2));
P = 1E-15; %photosynthetic rate in mol/cell/s
Pv = P ./ ((4*pi/3) * (par.Rc - par.st(1))^3); %volumentric photosynthetic rate mol/L/s NOTE that photosyntheis is negeletcted in the membrane layer so volume used is only the intracellular space
Pv = [Pv; Pv; Pv; 0; 0; 0; 0]; %put in vector form


%formulate reaction matrix in DBL (outside cell)
kuf = par.kp1 + par.kp4 .* par.OH;          %uncatalyzed CO2 hydration rate constant (/s)
kur = kuf .* par.bicfrac .* par.H./par.K1;  %uncatalzyed HCO3- dehydration rate (/s), removes fraction of "HCO3-" pool that is actually CO32-  
Ro = diag([-kuf, -kuf, -kuf, -kur, -kur, -kur, -kur]);
Ro(1,4) = kur; Ro(1,5) = (1/3)*kur; Ro(2,5) = (2/3) * kur; Ro(2,6) = (2/3)*kur; Ro(3,6) = (1/3) * kur; Ro(3,7) = kur;
Ro(4,1) = kuf; Ro(5,2) = kuf; Ro(6,3) = kuf;

%formulate reaction matrix inside cell
CAx = 1;        %factor by which CA is enhancing hyd/dehyd rates
Ri = CAx .* Ro;

%eCA activity
eCA = 1;        %factor by which eCA activity in surface layer enhances reaction rates.

%formulate diffusion matrix
dv = (4*pi) .* [par.Dco2, par.Dco2, par.Dco2, par.Db, par.Db, par.Db, par.Db]; %diffusivity vector, include 4*pi term
D = diag(dv);

%diffusion matrix for membrane layer
Dmem = diag(4*pi .* [par.Dco2, par.Dco2, par.Dco2, 0 , 0, 0, 0]);

%compute derivative vector in cell
for i = 1:par.shells(1)  
    
    R = Ri;             %chose correct reaction matrix
    j = 7*(i-1)+1;      %starting index
    k = 7*(i-1)+7;      %ending index, 7 species are treated (3 CO2, 4 HCO3-)

    %inner radius and outer radius of shell
    ri = par.st(1) .* (i-1);
    ro = par.st(1) .* i;
    
    %volume of shell (or sphere for inner-most layer)
    Vsh = (4*pi/3).*(ro^3 - ri^3);
    
    %distribute photosynthetic rate among 13CO2 species 
    ciFrac  = y(j:j+2,1)  ./ ([1 1 1] * y(j:j+2,1));
    Pvf = [ciFrac; 0; 0; 0; 0] .* Pv;
    
     % compute derivative vector (d[C]/dt, M/s)
    if (i == 1)
         ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).* (D *    (ro^2 .* (y(j+7:k+7,1) - y(j:k,1))./par.st(1))) - Pvf;  %inner-most sphere of cell
    elseif (i == par.shells(1)-1)  
         ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).* (Dmem * (ro^2 .* (y(j+7:k+7,1) - y(j:k,1))./par.st(1)) + D * (ri^2 .* (y(j-7:k-7,1) - y(j:k,1))./par.st(1))) - Pvf;   %shell just within membrane - account of lack of HCO3-/CO32- diffusion into membrane
    elseif (i == par.shells(1))
         ydot(j:k,1) =                  (1/Vsh).* (Dmem * (ro^2 .* (y(j+7:k+7,1) - y(j:k,1))./par.stInt       + ri^2 .* (y(j-7:k-7,1) - y(j:k,1))./par.st(1)));      %membrane layer, %treat as the cell membrane, permeable to CO2, but not to HCO3- or CO32-, neglect rxns in the membrane b/c rates and Ci concentrations (other than CO2)are ill-defined
    else
         ydot(j:k,1) = (R * y(j:k,1)) + (1/Vsh).* (D *    (ro^2 .* (y(j+7:k+7,1) - y(j:k,1))./par.st(1)       + ri^2 .* (y(j-7:k-7,1) - y(j:k,1))./par.st(1))) - Pvf;      %intermediate shells
    end
end

%compute derivative vector in DBL

for i = par.shells(1)+1:tshells
    
    R = Ro;             %chose correct reaction matrix
    j = 7*(i-1)+1;      %starting index
    k = 7*(i-1)+7;      %ending index, 3 species are treated (CO2, HCO3-, CO32-)
   
    %inner radius and outer radius of shell
    ri = par.Rc + par.st(2) .* (i - par.shells(1) - 1);
    ro = par.Rc + par.st(2) .* (i - par.shells(1));
    
    %volume of shell
    Vsh = (4*pi/3).*(ro^3 - ri^3);

     % compute derivative vector (d[C]/dt, M/s)
    if (i == par.shells(1)+1)
        ydot(j:k,1) = ((eCA .*R) * y(j:k,1)) + (1/Vsh).*(D * (ro^2 .* (y(j+7:k+7,1)     - y(j:k,1))./par.st(2)) + Dmem * (ri^2 .* (y(j-7:k-7,1) - y(j:k,1))./par.stInt));  % shell against cell surface
    elseif (i == tshells)
        ydot(j:k,1) = (       R  * y(j:k,1)) + (1/Vsh).*(D * (ro^2 .* (par.Cinit(1:7,1) - y(j:k,1))./par.st(2)          + ri^2 .* (y(j-7:k-7,1) - y(j:k,1))./par.st(2)));      %shell against bulk fluid
    else
        ydot(j:k,1) = (       R  * y(j:k,1)) + (1/Vsh).*(D * (ro^2 .* (y(j+7:k+7,1)     - y(j:k,1))./par.st(2)          + ri^2 .* (y(j-7:k-7,1) - y(j:k,1))./par.st(2)));      %intermediate shells
    end
end
    
ydot;
    
return


