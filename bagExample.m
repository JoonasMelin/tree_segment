%The path to your bag
bagPath = '/media/storage/20141007_evo_area1/2014-10-07-13-35-12_processed.bag';
bag = ros.Bag.load(bagPath);
bag.info()

topics = {'/copter2/camera/left/image_raw', ...
    '/copter2/camera/right/image_raw'};
startTime = 1412678168.44; %Change to the correct timestamp
bag.resetView(topics, startTime, []);

%Getting the tf tree
tree = ros.TFTree(bag);
tree.allFrames()

%Testing the tf
poseCamera = tree.lookup('earth_fixed', 'imu', startTime+startOffset);

%Testing the image rectification
bag.resetView(topics, startTime, []);
[msg, meta] = bag.read();

if strcmp(meta.topic, '/copter2/camera/right/image_raw') && ...
    timeRight == 0
    img2 = reshape(double(msg.data), [msg.width, msg.height])';
    img2_norm = mat2gray(img2);
    timeRight = (meta.time.time);
elseif strcmp(meta.topic, '/copter2/camera/left/image_raw') && ...
    timeLeft == 0
    img1 = reshape(double(msg.data), [msg.width, msg.height])';
    img1_norm = mat2gray(img1);
    timeLeft = (meta.time.time);
end
%Make sure to test that the images are taken at the same time!

%Rectifying the images
yErr = 2;
calibOk = 1;
while abs(yErr) > 1
    disp('Running online calibration..');
    try
        [R_cum, depthImg, rMin, rMax, img1_rect, img2_rect, dImg, yErr] = ...
            onlineStereoCalib(img1, img2, R_cum);
    catch
        disp('cannot calibrate');
        calibOk = 0;
        break;
    end
    disp(['Error left: ', num2str(yErr)]);
end