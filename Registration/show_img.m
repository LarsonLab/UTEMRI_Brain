function my_img = show_img(direc)
%Function shows the 3D image slices given a directory of dicom files.
dicomlist = dir(fullfile(direc,'*.DCM'));

files = cell(length(dicomlist),1);

for i = 1:length(files)
    files{i,1} = dicomlist(i).name;
end
files = natsortfiles(files);

cd(direc);
my_img = [];
for i = 1:numel(dicomlist)
    slice = dicomread(files{i});
    my_img = cat(3,my_img,slice);  

end 
imshow3D(my_img);

    
