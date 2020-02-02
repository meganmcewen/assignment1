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

%generate empty arrays
Px = zeros(n,1);
Py = zeros(n,1);
Vx = zeros(n,1);
Vy = zeros(n,1);
V_X = zeros(n,1);
V_Y = zeros(n,1);

%initialize particles
Px(:,1) = xrange*rand(n,1);
Py(:,1) = yrange*rand(n,1);

%calculate vTH
vTH = sqrt(2*k*T/m)

%generate random velocity
 randAngle = 2*pi*rand(n,1);
 Vx(:,1) = vTH * cos(randAngle);
 Vy(:,1) = vTH  *sin(randAngle);

%get Gaussian distribution for velocity components (scattering)
V_X(:,1) = vTH.*randn(n,1);
V_Y(:,1) = vTH.*randn(n,1);

%check average is close to vTH
Avg = sqrt(V_X.^2 + V_Y.^2);
VAvg = mean(Avg);


%plot histograms
figure(1)
hist(Avg,n)
title(['Distribution of Average Scattering Velocity. VAvg = ',num2str(VAvg),'.']);

figure(2)
hist(V_X,n);
title('Distribution of Scattering Velocities - Vx Component');

figure(3)
hist(V_Y,n);
title('Distribution of Scattering Velocities - Vy Component');

%time loop
timeStep = 1e-14;
time = 1:iter;

%set up temperature calculation
temp = zeros(iter,1);

%begin particle updating loop
for i = 1:iter
    
    %update dt
    dt = timeStep;

    %scattering
    Pscat = 1-exp(-dt/tau);
    ind = Pscat > rand(n,1);

    Vx(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
    Vy(ind) = sqrt((k*T)/m).*randn(sum(ind),1);

   %save old position
   Px_old = Px;
   Py_old = Py;
   
   %update position
   Px = Px + Vx*timeStep; %updates position in x
   Py = Py + Vy*timeStep; %updates position in y
    
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
   
   %temperature plot
   VAvg = mean(Vx.^2 + Vy.^2);
   temp(i) = (1/2)*(m*(VAvg))*(1/k);
   
   figure(4)
   plot(time,temp);
   title('Average Temperature in K');
   
   %mean free path
   MFP = VAvg*tau;
   
   %mean time between collisions
   MTC = (n./sum(VAvg))*MFP;
   
  %plot
  figure(5)
   plot([Px_old'; Px';] ,[Py_old'; Py';],"b");
   title(['Mean Free Path = ',num2str(MFP),' and Mean Time Between Collisions = ',num2str(MTC)]);
   axis ([0 200e-9 0 100e-9]);
   hold on
   drawnow

end
