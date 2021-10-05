function Flux = Cifluxes(Yss,p)
%calculate net diffusional fluxes and CO2/HCO3 interchange at steady state

ce = Yss(1,:);      %CO2 in extracellular solution
be = Yss(2,:);      %HCO3- in extracellular solution
cc = Yss(3,:);      %CO2 in cytoplasm
bc = Yss(4,:);      %HCO3- in cytoplasm
cp = Yss(5,:);      %CO2 in chloroplast stroma
bp = Yss(6,:);      %HCO3- in chloroplast stroma
cy = Yss(7,:);      %CO2 in pyrenoid
by = Yss(8,:);      %HCO3- in pyrenoid

%diffusional fluxes
Flux.cec = p.fc_c .*(ce - cc);      %CO2 flux from extracellular solution to cytoplasm
Flux.bec = p.fb_c .*(be - bc);      %HCO3- flux from extracellular solution to cytoplasm
Flux.ccp = p.fc_p .*(cc - cp);      %CO2 flux from cytoplasm into chloroplast stroma
Flux.bcp = p.fb_p .*(bc - bp);      %HCO3- flux from cytoplasm into chloroplast stroma
Flux.cpy = p.fc_y .*(cp - cy);      %CO2 flux from chloroplast stroma into pyrenoid
Flux.bpy = p.fb_y .*(bp - by);      %HCO3- flux from chloroplast stroma into pyrenoid

%CO2 dehydration fluxes
Flux.hyd_e    = (p.kuf .*ce).* (p.Ve./p.N);
Flux.dehyd_e  = (p.kur .*be).* (p.Ve./p.N);
Flux.hyd_c    = (p.kcf .*cc).* p.Vc;
Flux.dehyd_c  = (p.kcr .*bc).* p.Vc;
Flux.hyd_p    = (p.kpf .*cp).* p.Vp;
Flux.dehyd_p  = (p.kpr .*bp).* p.Vp;
Flux.hyd_y    = (p.kyf .*cy).* p.Vy;
Flux.dehyd_y  = (p.kyr .*by).* p.Vy;

Flux;
end