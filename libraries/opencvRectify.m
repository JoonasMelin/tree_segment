function rect = opencvRectify(img, K,K_new, R, d)
%The opencv implementation of the rectification algorithm
%http://docs.opencv.org/modules/imgproc/doc/geometric_transformations.html#void initUndistortRectifyMap(InputArray cameraMatrix, InputArray distCoeffs, InputArray R, InputArray newCameraMatrix, Size size, int m1type, OutputArray map1, OutputArray map2)
xSize = size(img,2);
ySize = size(img,1);
cx_hat  = K_new(1,3);
cy_hat = K_new(2,3);
fx_hat = K_new(1,1);
fy_hat = K_new(2,2);
fx = K(1,1);
fy = K(2,2);
cx = K(1,3);
cy = K(2,3);

[u,v] = meshgrid(1:xSize, 1:ySize);
x = (u-cx_hat)./fx_hat;
y = (v-cy_hat)./fy_hat;

cHomog = (inv(R)*[x(:), y(:), ones(length(x(:)), 1)]')';
cHomog = bsxfun(@rdivide,cHomog, cHomog(:,3)); %normalization
x_h = cHomog(:,1);
y_h = cHomog(:,2);

k1 = d(1);
k2 = d(2);
k3 = d(5);
p1 = d(3);
p2 = d(4);
r = sqrt((x_h.^2)+(y_h.^2));
one = ones(length(r),1);

x_hh = x_h.*(one+k1*r.^2+k2*r.^4+k3*r.^6) ...
    +2*p1*x_h.*y_h+p2*(r.^2+2*x_h.^2);

y_hh = y_h.*(one+k1*r.^2+k2*r.^4+k3*r.^6) ...
    +p1*(r.^2+2*y_h.^2)+2*(p2*x_h.*y_h);

map_x = bsxfun(@plus,(x_hh.*fx), cx);
map_y = bsxfun(@plus,(y_hh.*fy), cy);

%interpolating the new image
grid_x = reshape(map_x,size(img));
grid_y = reshape(map_y,size(img));
rect = interp2(u,v,double(img),grid_x, grid_y);
end