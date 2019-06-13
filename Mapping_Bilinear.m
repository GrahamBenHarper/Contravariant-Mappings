function point_physical = Mapping_Bilinear(cell_vertices,point_reference)

xhat = point_reference(1);
yhat = point_reference(2);

cell_x = cell_vertices(:,1);
cell_y = cell_vertices(:,2);

a(1) = cell_x(1);
a(2) = cell_x(2)-cell_x(1);
a(3) = cell_x(3)-cell_x(1);
a(4) = cell_x(4)-cell_x(3)-cell_x(2)+cell_x(1);

b(1) = cell_y(1);
b(2) = cell_y(2)-cell_y(1);
b(3) = cell_y(3)-cell_y(1);
b(4) = cell_y(4)-cell_y(3)-cell_y(2)+cell_y(1);

x = a(1) + a(2)*xhat + a(3)*yhat + a(4)*xhat*yhat;
y = b(1) + b(2)*xhat + b(3)*yhat + b(4)*xhat*yhat;

point_physical = [x,y];
return