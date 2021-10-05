function C = Phyto_1D(infile)
%simulate effects of photosynthesis and respiration on phtyoplankton surface C
%system and pH
%1-D spatial resolution (spherical system) from phytoplankton surface to bulk fluid

%load most parameters for the model
par = load1Dpar();

%initial system values at equilibrium which also serve as boundary values
%(i.e. bulk solution)
Cinit = repmat(par.Cinit,par.tshells,1);

%use ode solver to deterimine time course of isotope exchange
time = [0 100];
options = odeset('RelTol', 1E-6, 'AbsTol', 1E-9);         %set ode options to obtain smooth (non-oscillating or diverging) solution
[t_ode, Y] = ode15s(@D1Isoxderiv, time, Cinit, options, par);
Yt = Y';                  %ode gives one row per time point, transpose to return to normal structure: each row is the time series of a single Cspecies

%sample plot of species vs time in innermost layer (inside cell)
figure(1)
subplot(1,2,1)
plot(t_ode, Yt(1,:),'bo',t_ode, Yt(2,:),'go',t_ode, Yt(3,:),'ro'),title('CO2');
subplot(1,2,2)
plot(t_ode, Yt(4,:),'bo',t_ode, Yt(5,:),'go',t_ode, Yt(6,:),'ro',t_ode, Yt(7,:),'mo'),title('HCO3-');
[ax1,h3]=suplabel('DIC center of cell vs. Time','t');
set(h3,'FontSize',14) ;

figure(2)
subplot(1,2,1)
plot(t_ode, Yt(end-6,:),'bo',t_ode, Yt(end-5,:),'go',t_ode, Yt(end-4,:),'ro'),title('CO2 near bulk vs time');
subplot(1,2,2)
plot(t_ode, Yt(end-3,:),'bo',t_ode, Yt(end-2,:),'go',t_ode, Yt(end-1,:),'ro',t_ode, Yt(end,:),'mo'),title('HCO3-');
[ax1,h3]=suplabel('DIC near bulk of cell vs. Time','t');
set(h3,'FontSize',14) ;

%plot C species as a function of space at final timepoint - which should be
%the steady state solution
for k = 1:7
    Yf(k,:) = Y(end,k:7:end);
end
Yf = [Yf par.Cinit];

x1 = [par.st(1)/2:par.st(1):(par.Rc - par.st(1)/2)];        %from inner-most spherical compartment out to the membrane layer
x2 = [par.Rc:par.st(2):par.L];        %from surface DBL layer out to bulk solution, the xs are a bit off, but considered ok for now
x = [x1 x2];

%plot CO2 species as a function of distance from center of the cell out to
%bulk solution
figure(3)
subplot(1,2,1)
plot(x,Yf(1,:),'bo', x, Yf(2,:),'go', x, Yf(3,:),'ro'),title('CO2');
subplot(1,2,2)
plot(x,Yf(4,:),'bo', x, Yf(5,:),'go', x, Yf(6,:),'ro', x, Yf(7,:),'mo'),title('HCO3-');
[ax1,h3]=suplabel('DIC vs. x','t');
set(h3,'FontSize',14) ;


return

