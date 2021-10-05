function C = Phyto_1D(infile)
%simulate effects of photosynthesis and respiration on phtyoplankton surface C
%system and pH
%1-D spatial resolution (spherical system) from phytoplankton surface to bulk fluid

%load most parameters for the model
par = load1Dpar();

%initial system values at equilibrium which also serve as boundary values
%(i.e. bulk solution)
Cinit = [par.CO2; par.B; par.C];
par.Cinit = Cinit;
Cinit = repmat(Cinit,(par.shells(1) + par.shells(2)),1);

%use ode solver to deterimine time course of isotope exchange
time = [0 100];
options = odeset('RelTol', 1E-6, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Y] = ode15s(@D1deriv, time, Cinit, options, par);
Yt = Y';                  %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies

%sample plot of species vs time in innermost layer (inside cell)
figure(1)
subplot(1,3,1)
plot(t_ode, Yt(1,:)),title('CO2');
subplot(1,3,2)
plot(t_ode, Yt(2,:)),title('HCO3-');
subplot(1,3,3)
plot(t_ode, Yt(3,:)), title('CO32-');

%plot C species as a function of space at final timepoint - which should be
%the steady state solution
Yf(1,:) = Y(end,1:3:end);           %CO2
Yf(2,:) = Y(end,2:3:end);           %HCO3-
Yf(3,:) = Y(end,3:3:end);           %CO32-
Yf = [Yf Cinit(1:3,1)];

x1 = [par.st(1)/2:par.st(1):(par.Rc - par.st(1)/2)];        %from inner-most spherical compartment out to the membrane layer
x2 = [par.Rc:par.st(2):par.L];        %from surface DBL layer out to bulk solution, the xs are a bit off, but considered ok for now
x = [x1 x2];

%plot CO2 species as a function of distance from center of the cell out to
%bulk solution
figure(2)
subplot(1,3,1)
plot(x,Yf(1,:),'ro'),title('CO2');
subplot(1,3,2)
plot(x,Yf(2,:),'go'),title('HCO3-');
subplot(1,3,3)
plot(x,Yf(3,:),'bo'),title('CO32-');

%plot just the data inside the cell

figure(3)
subplot(1,3,1)
plot(x1,Yf(1,1:par.shells(1)),'ro'),title('CO2 inside cell');
subplot(1,3,2)
plot(x1,Yf(2,1:par.shells(1)),'go'),title('HCO3-');
subplot(1,3,3)
plot(x1,Yf(3,1:par.shells(1)),'bo'),title('CO32-');


return

