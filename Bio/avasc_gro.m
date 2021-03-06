%avascular tumor growth model
xsize = 100; %size of matrix
ysize = 100;
temp = zeros(xsize,ysize);  %for storing and assigning updated values
c = zeros(xsize,ysize); % concentration of nutrients
c0 = rand;   %the original concentration of nutrients c naught, randomly generated between 0 & 1
c(1,1) = c0;
xmax = 1;    %x coord of boundary 
ymax = 1;    %y coord of boundary 
n_iter = 250;    %number of iterations
dx = xmax/(xsize-1);   %array spacing
dy = ymax/(ysize-1);
x = 0:dx:xmax;
y = 0:dy:ymax;
f = zeros(xsize,ysize); %function of nutrient concentration
g = zeros(xsize,ysize); %func of c
h = zeros(xsize,ysize); %func of c
p = zeros(xsize,ysize); %proliferating cells
q = zeros(xsize,ysize); %quiscient cells
n = zeros(xsize,ysize); %necrotic cells
r = zeros(xsize,ysize); %living cells
u = zeros(xsize,ysize); %func of c
v = zeros(xsize,ysize); %func of c
p(:,1) = 1; %implement boundary conditions
p(1,:) = 1;
q(:,1) = 0;
n(:,1) = 0;
c(1,:) = c0;
c(:,1) = 1;
alph = 0.9; %constants
beta = 0.5;
gamm = 10;
t = linspace(0,10); %time
dt = 10/99;
diff = rand;
i = 2:xsize-1;    %define indexing integers for iteration
j = 2:ysize-1;
for it=1:n_iter
    f(i,j) = 0.5*(1-tanh(4*c(i,j)-2));
    h(i,j) = 0.25*(1-tanh(4*c(i,j)-2));
    g(i,j) = 0.5*exp(0.5*c(i,j));
    p(i,j) = p(i,j)+dt*(u(i,j)+g(i,j)*p(i,j)*(1-p(i,j)-q(i,j)-n(i,j))-f(i,j)*p(i-1,j-1));
    q(i,j) = q(i,j)+dt*(v(i,j)+f(i,j)*q(i,j)-h(i,j)*q(i,j));
    n(i,j) = n(i,j)+dt*(h(i,j)*q(i,j));
    c(i,j) = (gamm./(gamm+p(i,j)))*(1-alph.*(p(i,j)+q(i,j)+n(i,j)));
    r(i,j) = p(i-1,j-1)+q(i-1,j-1);
    u(i,j) = ((p(i+1,j)-p(i-1,j))*r(i,j)*(r(i+1,j)-r(i-1,j))+4*p(i,j)*r(i,j)*(r(i+1,j)-2*r(i,j)+r(i-1,j))-p(i,j)*(r(i+1,j)-r(i-1,j))^2)/(4*dx^2*r(i,j)^2);
    v(i,j) = ((q(i+1,j)-q(i-1,j))*r(i,j)*(r(i+1,j)-r(i-1,j))+4*q(i,j)*r(i,j)*(r(i+1,j)-2*r(i,j)+r(i-1,j))-q(i,j)*(r(i+1,j)-r(i-1,j))^2)/(4*dx^2*r(i,j)^2);
    r(i,j) = p(i,j)+q(i,j);
end
%plot(t,p)

%Data needed to calibrate program, plots cell growth rates. central
%difference and linear difference approximations