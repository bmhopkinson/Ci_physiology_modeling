function ydot = Cideriv(t,y,p)
%derivative function for Ci species in mechanistic model of the CCM

ce = y(1,1);
be = y(2,1);
cc = y(3,1);
bc = y(4,1);
cp = y(5,1);
bp = y(6,1);
cy = y(7,1);
by = y(8,1);

%calculate uptake and photosynthetic rates
ind = find(p.ON <=t,1,'last');  %determine if light is on or off
if( isempty(ind) || t > p.OFF(ind) )          %light  off
    Bupc = 0;
    Bupp = 0;
    P    = 0;
else                        %light on
    Bupc  = (p.Vm_Bc.*be)./(p.Km_Bc + be);
    Bupp  = (p.Vm_Bp.*bc)./(p.Km_Bp + bc);
    P     = p.mRub.*(p.kcat_R .* cy)./(p.Km_R + cy);

end

dce = -p.kuf.*ce  + p.kur.*be  + (p.N./p.Ve).*p.fc_c.*(cc - ce);
dbe =  p.kuf.*ce  - p.kur.*be  + (p.N./p.Ve).*(p.fb_c.*(bc - be) - Bupc);
dcc = -p.kcf.*cc  + p.kcr.*bc  + (1/p.Vc).*(p.fc_c.*(ce - cc) + p.fc_p.*(cp - cc));
dbc =  p.kcf.*cc  - p.kcr.*bc  + (1/p.Vc).*(p.fb_c.*(be - bc) + p.fb_p.*(bp - bc) + Bupc - Bupp);
dcp = -p.kpf.*cp  + p.kpr.*bp  + (1/p.Vp).*(p.fc_p.*(cc - cp) + p.fc_y.*(cy - cp));
dbp =  p.kpf.*cp  - p.kpr.*bp  + (1/p.Vp).*(p.fb_p.*(bc - bp) + p.fb_y.*(by - bp) + Bupp);
dcy = -p.kyf.*cy  + p.kyr.*by  + (1/p.Vy).*(p.fc_y.*(cp - cy) - P);
dby =  p.kyf.*cy  - p.kyr.*by  + (1/p.Vy).*(p.fb_y.*(bp - by));

ydot = [dce; dbe; dcc; dbc; dcp; dbp; dcy; dby];

end

