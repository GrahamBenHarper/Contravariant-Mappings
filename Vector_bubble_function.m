function phi = Vector_bubble_function(point)

xhat = point(1);
yhat = point(2);

normal = [1,0];

phi = 4*normal*xhat*yhat*(1-yhat);
return