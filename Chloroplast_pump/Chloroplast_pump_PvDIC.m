function out = Chloroplast_pump_PvDIC(infile)
%mechanistic model of the CCM to match P vs DIC data

par1file = strcat(infile,'.par');
par2file = strcat(infile,'.par2');
CO2dat  = strcat(infile,'_cdat.txt');   %calibrated CO2 data from PvDIC experiment
[p, Fobs] = loadparams1(par1file);     %experimental parameters and rates derived from O2 and CO2 data
p = loadparams2(p, par2file);          %model parameters

p.DICback = 25E-9;      %background 13CDIC concentration to for initial conditions
yinit = initcond(p);

fdat = fopen(CO2dat,'r');
C = textscan(fdat,'%f %f');
tdat = C{1};
CO2dat = C{2};

tbrks =[0 p.ADD p.OFF(end-1)];

t = [];
Y = [];
%use ode solver to deterimine time course of CO2 and HCO3- behavior
for i  = 1:length(p.ADD)+1
    tspan = tbrks(i:i+1);
    options = odeset('RelTol', 1E-6, 'AbsTol', 1E-10,'MaxStep',5);         %set ode options to obtain smooth (non-oscillating or diverging) solution
    [t_ode, Ys] = ode15s(@Cideriv, tspan, yinit, options, p);
    Ys = Ys';                  %ode gives one row per time point, transpose structure in which each row is the time series of a single Cspecies
    Y = [Y Ys];             %append new data onto final solution array
    t = [t t_ode'];
    if (i <= length(p.DIC))
        yinit = Ys(:,end);
        yinit(2,1) = yinit(2,1) + p.DIC(i);
    end

end
    
figure(1)    
subplot(2,4,1)
plot(t,Y(1,:),'r'), title('CO2_e');
subplot(2,4,5)
plot(t,Y(2,:),'r'), title('HCO3_e');
subplot(2,4,2)
plot(t,Y(3,:),'r'),title('CO2_c');
subplot(2,4,6)
plot(t,Y(4,:),'b'), title('HCO3_c');
subplot(2,4,3)
plot(t,Y(5,:),'b'),title('CO2_p');
subplot(2,4,7)
plot(t,Y(6,:),'b'), title('HCO3_p');
subplot(2,4,4)
plot(t,Y(7,:),'b'),title('CO2_y');
subplot(2,4,8)
plot(t,Y(8,:),'b'), title('HCO3_y');

figure(2)
plot(t,Y(1,:),'r',tdat,CO2dat,'og'), title('CO2_e');

%compute determined P, Uptake rates
ind2 = arrayfun(@(x) find(t < x,1,'last'), p.OFF(1:end-2));
stop = find(t < p.OFF(end-1),1,'last');
ind2 = [ind2 stop];

Yss = zeros(8,length(ind2));
for i = 1:length(ind2)
    Yss(:,i) = mean(Y(:,ind2(i)-6:ind2(i)-1),2);
end

P = p.mRub.*(p.kcat_R .* Yss(7,:))./(p.Km_R + Yss(7,:));        %photosynthetic rate
Cup_c = p.fc_c .* (Yss(1,:) - Yss(3,:));        %net CO2 uptake
Bup_c = (p.Vm_Bc .* Yss(2,:))./(p.Km_Bc + Yss(2,:));    %HCO3- uptake into cell
Bup_p = (p.Vm_Bp .* Yss(4,:))./(p.Km_Bp + Yss(4,:));    %chloroplast pump HCO3- transport rate

%calculate fluxes 
Flux = Cifluxes(Yss,p);

figure(3)
plot(Yss(2,:),P,'r',Yss(2,:),Cup_c,'g',Yss(2,:),p.kur.*Yss(2,:)./p.N,'-.g',Yss(2,:), Bup_c,'b',Yss(2,:), Bup_p,'m',Yss(2,:),Fobs.P,'or',Yss(2,:),Fobs.Cup,'og', Yss(2,:), Fobs.Bup,'ob');

%write data out to files
dlmwrite('Chloroplast_pump_Fulloutput.txt',[t' Y'],'\t');           %full time course of Ci concentrations in all compartments

ssfile = 'Chloroplast_pump_Fluxes.txt';
fid = fopen(ssfile,'w');
fprintf(fid,'CO2e\t HCO3e\t CO2c\t HCO3c\t CO2p\t HCO3p\t CO2y\t HCO3y\t P\t Bup_c\t Bup_y\t Diff_CO2_etoc\t Diff_HCO3_etoc\t Diff_CO2_ctop\t Diff_HCO3ctop\t Diff_CO2_ptoy\t Diff_HCO3_ptoy\t Hyd_e\t Dehyd_e\t Hyd_c\t Dehyd_c\t Hyd_p\t Dehyd_p\t Hyd_y\t Dehyd_y\n'); 
fclose(fid);
dlmwrite(ssfile,[Yss' P' Bup_c' Bup_p' Flux.cec' Flux.bec' Flux.ccp' Flux.bcp' Flux.cpy' Flux.bpy' Flux.hyd_e' Flux.dehyd_e' Flux.hyd_c' Flux.dehyd_c' Flux.hyd_p' Flux.dehyd_p' Flux.hyd_y' Flux.dehyd_y'],'-append','delimiter','\t');

end

