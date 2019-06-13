function J = Jacobian_Bilinear(cell_vertices,point_reference)

J = zeros(2,2);

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

J(1,1) = a(2) + a(4)*yhat;
J(1,2) = a(3) + a(4)*xhat;
J(2,1) = b(2) + b(4)*yhat;
J(2,2) = b(3) + b(4)*xhat;

return