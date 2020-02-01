%constants
xrange = 2e-7; %size of area in x
yrange = 1e-7; %size of area in y
n = 10; %number of particles
m0 = 9.10938356e-31; %electron mass
m = 0.26*m0;
T = 300; %temperature (K)
k = 1.380648e-23; %Boltzmann constant
tau = 0.2e-12;
iter = 500; %number of iterations to run the simulation

%calculate vTH
vTH = sqrt(2*k*T/m);

%generate empty arrays
Px = zeros(n,1);
Py = zeros(n,1);
Vx = zeros(n,1);
Vy = zeros(n,1);

%initialize particles
Px(:,1) = xrange*rand(n,1);
Py(:,1) = yrange*rand(n,1);

%generate random velocity
randAngle = 2*pi*rand(n,1);
Vx(:,1) = vTH * cos(randAngle);
Vy(:,1) = vTH  *sin(randAngle);

%time loop
timeStep = 1e-14;
time = 1:iter;

%set up temperature calculation
temp = zeros(iter,1);

%begin particle updating loop
for i = 1 : iter
    
   %save old position
   Px_old = Px;
   Py_old = Py;
   
   %update position
   Px = Px + Vx*timeStep; %updates position in x
   Py = Py + Vy*timeStep; %updates position in y
    
  %calculate mean free path
  VAvg = mean(Vx.^2 + Vy.^2);
  MFP = VAvg*tau
  
  %calculate temperature
  temp(i) = (1/2)*(m*(VAvg))*(1/k);
  
   %boundaries
   %x hitting right side
   id = Px >= xrange;
   Px(id) = Px(id) - xrange;
   Px_old(id) = Px_old(id) - xrange;

   %x hitting left side
   id = Px <= 0;
   Px(id) = Px(id) + xrange;
   Px_old(id) = Px_old(id) + xrange;
 
   %bouncing y off top/bottom
   Vy(Py >= yrange) = Vy(Py >= yrange) * -1;
   Vy(Py <= 0) = Vy(Py <= 0) * -1;
   Py(Py>yrange) = yrange-(Py(Py>yrange)-yrange);

  %plot particles
  figure(1)
   plot([Px_old'; Px';] ,[Py_old'; Py';],"b");
   axis ([0 200e-9 0 100e-9]);
   hold on
   drawnow
   
   %plot temperature
   figure(2)
   plot(time,temp);
   title('Average Temperature in K');
   
end


