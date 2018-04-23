xsize = 100; %size of matrix
ysize = 100;
temp = zeros(xsize,ysize);  %for storing and assigning updated values
c = zeros(xsize,ysize); %final potential solution
c0 = rand;   %the concentration of VEGF c naught, randomly generated between 0 & 1
xmax = 1;    %x coord of boundary (width of shape/area)
ymax = 1;    %y coord of boundary (height of shape/area)
n_iter = 250;    %number of iterations
dx = xmax/(xsize-1);   %array spacing
dy = ymax/(ysize-1);
x = 0:dx:xmax;
y = 0:dy:ymax;
temp(:,1) = 0; %implement boundary conditions (increasing chem gradient)
temp(:,ysize) = 0;
temp(1,:) = c0;
temp(xsize,:) = temp(xsize-1,:);
i = 2:xsize-1;    %define indexing integers for iteration
j = 2:ysize-1;
for it=1:n_iter
    c = temp;
    temp(i,j) = (c(i+1,j) + c(i-1,j) + c(i,j+1) + c(i,j-1))/4;  %average neighbors
    temp(:,1) = 0; %re-implement boundary conditions
    temp(:,ysize) = 0;
    temp(1,:) = c0;
    temp(xsize,:) = temp(xsize-1,:);
end
t = linspace(0,10);
dt = 10/99;
diff = rand;
%c = [x.*rand; y.*rand];
p0 = zeros(100,2);
p1 = zeros(100,2);
p2 = zeros(100,2);
p3 = zeros(100,2);
p4 = zeros(100,2);
v1 = zeros(100,2);
v2 = zeros(100,2);
v3 = zeros(100,2);
v4 = zeros(100,2);
v2(:,1) = 0.25;
v3(:,1) = 0.5;
v4(:,1) = 0.95;
chemFact = 10.*c;
%cgrad = gradient(c,[x,y]);
i = 2:99;
j = 2:99;
for it=1:100
    p1(i,j) = (dt*diff/(dx^2))-(dt/(4*dx^2))*(chemFact(i,j)*(c(i+1,j)-c(i-1,j)));
    p2(i,j) = (dt*diff/(dx^2))+(dt/(4*dx^2))*(chemFact(i,j)*(c(i+1,j)-c(i-1,j)));
    p3(i,j) = dt*diff/(dx^2)-(dt/(4*dx^2))*(chemFact(i,j)*(c(i,j+1)-c(i,j-1)));
    p4(i,j) = dt*diff/(dx^2)+(dt/(4*dx^2))*(chemFact(i,j)*(c(i,j+1)-c(i,j-1)));
    p0(i,j) = 1 - (p1(i,j) + p2(i,j) + p3(i,j) + p4(i,j));
    p0(:,100) = 1;
    p1(:,100) = 0;
    p2(:,100) = 0;
    p3(:,100) = 0;
    p4(:,100) = 0;
end
for it=1:n_iter
    if p1(i,j)>p0(i,j)
        v1(i+1,1) = v1(i,1) - dx;
    end
    if p2(i,j)>p0(i,j)
        v1(i+1,1) = v1(i,1) + dx;
    end
    if p3(i,j)>p0(i,j)
        v1(j+1,2) = v1(j,2) - dy;
    end
    if p4(i,j)>p0(i,j)
        v1(j+1,2) = v1(j,2) + dy;
    else
        v1(i,1) = v1(i,1);
        v1(i,2) = v1(i,2);
    end
end
v1x = v1(:,1);
v1y = v1(:,2);
