%constants
xrange = 2e-7; %size of area in x
yrange = 1e-7; %size of area in y
n = 100; %number of particles
m0 = 9.10938356e-31; %electron mass
m = 0.26*m0;
T = 300; %temperature (K)
k = 1.380648e-23; %Boltzmann constant
tau = 0.2e-12;
iter = 100; %number of iterations to run the simulation

%generate empty arrays for particles
Px = zeros(n,1);
Py = zeros(n,1);
Vx = zeros(n,1);
Vy = zeros(n,1);
V_X = zeros(n,1);
V_Y = zeros(n,1);

%initialize particle positions
Px(:,1) = xrange*rand(n,1);
Py(:,1) = yrange*rand(n,1);

%calculate vTH
vTH = sqrt(2*k*T/m);

%generate random velocity
randAngle = 2*pi*rand(n,1);
Vx(:,1) = vTH * cos(randAngle);
Vy(:,1) = vTH  *sin(randAngle);

%get Gaussian distribution for velocity components (scattering)
V_X(:,1) = vTH.*randn(n,1);
V_Y(:,1) = vTH.*randn(n,1);

%check average of distribution is close to vTH
Avg = sqrt(V_X.^2 + V_Y.^2);

%plot histograms
figure(1)
hist(Avg,n);
title('Distribution of Average Scattering Velocity');

%time loop setup
timeStep = 1e-14;
time = 1:iter;

%set up temperature calculation
temp = zeros(iter,1);

%define boxes for bottleneck
boxLeft = 0.8e-7;
boxRight = 1.2e-7;
boxBottom = 0.6e-7;
boxTop = 0.4e-7;

%define inside of boxes
inX = (boxLeft < Px) & (Px < boxRight);
inY = (boxBottom < Py) | (boxTop > Py);
inBox = inX&inY;

%check particles are not initially in boxes and move if needed
numInBox = 1;
while(numInBox > 0)
   
    Px(inBox) = xrange*rand(1,numInBox); 
  
    inX = (boxLeft < Px) & (Px < boxRight);
    inY = boxBottom < Py | boxTop > Py;
    inBox = inX&inY;
    
    numInBox = sum(inBox)
end
  
%begin particle updating loop
for i = 1 : iter
    
    %set dt
    dt = timeStep;
   
    %scattering
    Pscat = 1-exp(-dt/tau);
    
    ind = Pscat > rand(n,1);
    Vx(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
    Vy(ind) = sqrt((k*T)/m).*randn(sum(ind),1);
    
   %save old position of particles
   Px_old = Px;
   Py_old = Py;
   
   %update position
   Px = Px + Vx*timeStep; %updates position in x
   Py = Py + Vy*timeStep; %updates position in y
   
   %bouncing off boxes if box boundary encountered
   inX = (boxLeft < Px) & (Px < boxRight);
   inY = boxBottom < Py | boxTop > Py;
   inBox = inX&inY;
    
   betweenBoxes = (Py_old > boxTop)&(Py_old < boxBottom);
   Vy(inBox&betweenBoxes) = -1*Vy(inBox&betweenBoxes);
   Vx(inBox&~betweenBoxes) = -1*Vx(inBox&~betweenBoxes);
   
   %boundaries of total area
   %x hitting right side
   id = Px >= xrange;
   Px(id) = Px(id) - xrange;
   Px_old(id) = Px_old(id) - xrange;

   %x hitting left side
   id = Px <= 0;
   Px(id) = Px(id) + xrange;
   Px_old(id) = Px_old(id) + xrange;
 
   %bouncing y off top/bottom
   id = Py >= yrange;
   Vy(id) = Vy(id) * -1;
   Vy(Py <= 0) = Vy(Py <= 0) * -1;
   Py(id) = yrange-(Py(id)-yrange);
  
   %temperature plot
   VAvg = mean(Vx.^2 + Vy.^2);
   temp(i) = (1/2)*(m*(VAvg))*(1/k);
   
   figure(3)
   plot(time,temp);
   title('Average Temperature in K');
   
   %mean free path
   MFP = VAvg*tau;
   
  %plot particles
   figure(1)
   plot([Px_old'; Px';] ,[Py_old'; Py';],'b');
   axis ([0 200e-9 0 100e-9]);
   title('Particle Movement in a Defined Space with Scattering');
   hold on
   drawnow
   pause(0.05)
   %add boxes to particle plot
   rectangle('Position',[0.8e-7 0 0.4e-7 0.4e-7])
   rectangle('Position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7])
   
end

%electron density map with final positions
figure(4)
hist3([Px Py],'CdataMode','auto')
title(['Electron Density Map after ', num2str(iter),' iterations for ' ,num2str(n),' particles']);
colorbar
view(2)