% Create a VideoReader object for the video file
fileName = "db0edda2_0.0-1164.705.mp4";
vidObj = VideoReader(fileName);


numFrames = floor(vidObj.Duration * vidObj.FrameRate);
%%
vidObj.CurrentTime = 0; % Start from the beginning
fIdx = 0;
v=[];
% Read and display frames until the end of the video
while hasFrame(vidObj)
    vidFrame = readFrame(vidObj);
    % cat(4,v,vidFrame);
    fIdx = fIdx+1;
    disp(fIdx)
end
%%

% Create a figure to display the video frames
figure;

vidObj.CurrentTime = 0; % Start from the beginning
% Read the next frame
vidFrame = readFrame(vidObj);
imshow(vidFrame);

