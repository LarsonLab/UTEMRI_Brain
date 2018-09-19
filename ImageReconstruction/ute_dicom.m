function ute_dicom(finalImage, pfile_name, output_image, image_option, scaleFactor)
% Convert matlab 3D matrix to dicom for UWUTE sequence
% resolution is fixed in the recon - FOV/readout(from scanner), isotropic
% matrix size is determined in the recon
% Inputs:
%   finalImage: 3D image matrix
%   pfile_name: original pfile name
%   output_dir: output directory
%   image_option: 1 for both phase and magnitude, 0(offset) mag only
%
% August, 2018, Xucheng Zhu


if nargin<4
    image_option = 0;
end
addpath(genpath('../util'));


addpath(genpath('../orchestra-sdk-1.7-1.matlab'));
Isize = size(finalImage);

pfile = GERecon('Pfile.Load', pfile_name);
pfile.phases = numel(finalImage(1,1,1,1,:)); 
pfile.xRes = size(finalImage,1);
pfile.yRes = size(finalImage,2);
pfile.slices = size(finalImage,3);
pfile.echoes = size(finalImage,4);

% calc real res(isotropic,axial)
corners = GERecon('Pfile.Corners', 1);
res = abs(corners.UpperRight(2)-corners.UpperLeft(2))/Isize(3); 
scale = Isize/Isize(3); 
corners.LowerLeft = corners.LowerLeft.*scale;
corners.UpperLeft = corners.UpperLeft.*scale;
corners.UpperRight = corners.UpperRight.*scale;
orientation = GERecon('Pfile.Orientation', 1);

%  X = repmat(int16(0), [96 86 1 94]);
X = zeros(96, 86, 1, 94);


seriesNumber = randi(99999);
seriesDescription = ['UTE T2 - ', output_image];

for s = 1:pfile.slices
    
    for e = 1:pfile.echoes
        for p = 1:pfile.phases
            
            mag_t =double(finalImage(:,:,s,e,p) * scaleFactor);
%             figure;imshow(mag_t);title('mag_t');

%             mag_t = GERecon('Orient', mag_t, orientation);

            imageNumber = ImageNumber(s, e, p, pfile);
            filename = ['DICOMS/',output_image, '/image_',num2str(imageNumber) '.dcm'];
            GERecon('Dicom.Write', filename, mag_t, imageNumber, orientation, corners, seriesNumber, seriesDescription);
            if image_option~=0
                phase_t = flip(flip(single(angle(finalImage(:,:,s,e,p))).',1),2);
                %phase_t = GERecon('Orient', phase_t, orientation);
                filename = [output_dir,'DICOMS/phase_',num2str(imageNumber) '.dcm'];
                GERecon('Dicom.Write', filename, phase_t, imageNumber, orientation, corners);
            end
            
            [X(:,:,1,s),map] = dicomread(filename);
        end
    end
    
    % Get corners and orientation for this slice location
    corners.LowerLeft(3) = corners.LowerLeft(3) + res;
    corners.UpperLeft(3) = corners.UpperLeft(3) + res;
    corners.UpperRight(3) = corners.UpperRight(3) + res;
end

disp([output_image, ' generated.']);

figure;montage(X(:,:,1,:),map);title(output_image);

end
               

function number = ImageNumber(slice, echo, phase, pfile)
% Image numbering scheme:
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn
    slicesPerPhase = pfile.slices * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1) + 1;
end
