%% Setup
clear all
clf

% physical cell properties
%cell_vertices = [1,0; 2,1; 0,1; 1,2];
cell_vertices = [0.1,0; 0.2,0.1; 0,0.1; 0.1,0.2];

nx = 20;
ny = 20;

x = linspace(0,1,nx+1);
y = linspace(0,1,ny+1);

[X,Y] = meshgrid(x,y);
U = zeros(size(X));
V = zeros(size(Y));

%% reference cell
phi = zeros((nx+1)*(ny+1),2);
for i=1:length(x)
  for j=1:length(y)
    phi = Vector_bubble_function([X(i,j),Y(i,j)]);
    U(i,j) = phi(1);
    V(i,j) = phi(2);
  end
end

disp('Length on reference cell: 1')

figure(1)
quiver(X,Y,U,V)
title('Reference cell function')
axis([0,1,0,1])

%% physical cell with no transform
[X,Y] = meshgrid(x,y);
U = zeros(size(X));
V = zeros(size(Y));
maxlength = 0;

phi = zeros((nx+1)*(ny+1),2);
for i=1:length(x)
  for j=1:length(y)
    phi = Vector_bubble_function([X(i,j),Y(i,j)]);
    
    physical_point = Mapping_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    X(i,j) = physical_point(1);
    Y(i,j) = physical_point(2);
    U(i,j) = phi(1);
    V(i,j) = phi(2);
    
    l = sqrt( phi(1)^2 + phi(2)^2 );
    if(l>maxlength)
      maxlength=l;
    end
  end
end

disp(['Length with no transform: ' num2str(maxlength)])

figure(2)
quiver(X,Y,U,V)
hold on
plot(cell_vertices([1,2,4,3,1],1),cell_vertices([1,2,4,3,1],2),'k')
title('Physical cell function')

%% physical cell via contravariant transform
[X,Y] = meshgrid(x,y);
U = zeros(size(X));
V = zeros(size(Y));
maxlength = 0;

phi = zeros((nx+1)*(ny+1),2);
for i=1:length(x)
  for j=1:length(y)
    phi = Vector_bubble_function([X(i,j),Y(i,j)]);
    J = Jacobian_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    
    phi = (J*phi')';
    
    physical_point = Mapping_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    X(i,j) = physical_point(1);
    Y(i,j) = physical_point(2);
    U(i,j) = phi(1);
    V(i,j) = phi(2);
    
    l = sqrt( phi(1)^2 + phi(2)^2 );
    if(l>maxlength)
      maxlength=l;
    end
  end
end

disp(['Length with contravariant transform: ' num2str(maxlength)])

figure(3)
quiver(X,Y,U,V)
hold on
plot(cell_vertices([1,2,4,3,1],1),cell_vertices([1,2,4,3,1],2),'k')
title('Physical cell function via contravariant transform')

%% physical cell via Piola transform
[X,Y] = meshgrid(x,y);
U = zeros(size(X));
V = zeros(size(Y));
maxlength = 0;

phi = zeros((nx+1)*(ny+1),2);
for i=1:length(x)
  for j=1:length(y)
    phi = Vector_bubble_function([X(i,j),Y(i,j)]);
    J = Jacobian_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    
    phi = 1/det(J)*(J*phi')';
    
    physical_point = Mapping_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    X(i,j) = physical_point(1);
    Y(i,j) = physical_point(2);
    U(i,j) = phi(1);
    V(i,j) = phi(2);
    
    l = sqrt( phi(1)^2 + phi(2)^2 );
    if(l>maxlength)
      maxlength=l;
    end
  end
end

disp(['Length with Piola transform: ' num2str(maxlength)])

figure(4)
quiver(X,Y,U,V)
hold on
plot(cell_vertices([1,2,4,3,1],1),cell_vertices([1,2,4,3,1],2),'k')
title('Physical cell function via Piola transform')

%% physical cell via orthogonal contravariant transform
[X,Y] = meshgrid(x,y);
U = zeros(size(X));
V = zeros(size(Y));
maxlength = 0;

phi = zeros((nx+1)*(ny+1),2);
for i=1:length(x)
  for j=1:length(y)
    phi = Vector_bubble_function([X(i,j),Y(i,j)]);
    J = Jacobian_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    
    phi = 1/sqrt(det(J))*(J*phi')';
    
    physical_point = Mapping_Bilinear(cell_vertices,[X(i,j),Y(i,j)]);
    X(i,j) = physical_point(1);
    Y(i,j) = physical_point(2);
    U(i,j) = phi(1);
    V(i,j) = phi(2);
    
    l = sqrt( phi(1)^2 + phi(2)^2 );
    if(l>maxlength)
      maxlength=l;
    end
  end
end

disp(['Length with orthogonal contravariant transform: ' num2str(maxlength)])

figure(5)
quiver(X,Y,U,V)
hold on
plot(cell_vertices([1,2,4,3,1],1),cell_vertices([1,2,4,3,1],2),'k')
title('Physical cell function via orthogonal contravariant transform')