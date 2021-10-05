function [p, Fobs] = loadparams1(pFile)
%load pameters for P vs DIC experiments

fp = fopen(pFile,'r');

if fp==-1
   error(['File ' pFile ' not found or permission denied.']);
end

i=1;

%define pameter structure; initialize values that should be filled from
%pameter file
p = struct('N',[],'O2calib',[],'O2back',[],'CO2calib',[],'CO2back',[],...
    'pH',[],'TC',[]);

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
        
        case 'cells/mL'
            p.N = sscanf(raw{i},'%*s %e');        %number of cells/mL
        case 'O2calib'
            p.O2calib = sscanf(raw{i},'%*s %e');      %mol O2/cm3/A
        case 'O2background'
            p.O2back = sscanf(raw{i},'%*s %e');       %O2 background current (A) 
        case 'CO2calib'
            p.CO2calib = sscanf(raw{i},'%*s %e');     %mol CO2/A
        case 'CO2background'
            p.CO2back = sscanf(raw{i}, '%*s %i %i');  %cycle numbers for CO2 background
        case 'pH'
            p.pH = sscanf(raw{i},'%*s %f');           %pH of assay buffer
        case 'temp'
            p.TC = sscanf(raw{i},'%*s %f');        %temperature (in C on file, now in K)
    end
end

%check to make sure all the fields were filled
fn = fieldnames(p);
for i = 1:length(fn)
    if isempty(p.(fn{i}))
        error('Error in load_params, missing pameter: %s\n',fn{i})
    end
end


%read in light on/off cycle numbers
line = fgetl(fp);     %discard header line
i = 1;
while 1
    line = fgetl(fp);
    if (strcmp('________',line))
        break
    end
    A = sscanf(line, '%f %f %f');
    p.DIC(i) = A(1)./1E9;     %convert from uM to mol/cm3
    p.ON(i) = A(2);
    p.OFF(i) = A(3);
    i = i+1;
end
p.DIC = p.DIC(1:end-1);     %remove -1 at end of DIC vector

%read in data on DIC addition times and time that the CO2 signal stabilizes
line = fgetl(fp);     %discard header line
i = 1;
while 1
    line = fgetl(fp);
    if (strcmp('________',line))
        break
    end
    A = sscanf(line, '%f %f %f');
    p.ADD(i) = A(2);
    p.STA(i) = A(3);
    i = i+1;
end

%read in observed rates of photosynthesis and CO2 uptake
line = fgetl(fp);     %discard header line
i = 1;
while 1
    line = fgetl(fp);
    if (strcmp('________',line))
        break
    end
    A = sscanf(line, '%f %e %e');
    Fobs.P(i)   = A(2);
    Fobs.Cup(i) = A(3);
    i = i+1;
end
Fobs.Bup = Fobs.P - Fobs.Cup;       %calculate net HCO3- uptake as difference between photosynthesis and CO2 uptake

fclose(fp);

return
    

