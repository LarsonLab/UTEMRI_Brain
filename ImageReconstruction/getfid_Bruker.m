function cplx = getfid(filenm, dim, rephase)

if nargin == 0
    error('No fid specified. Give at least 1 argument');
end

if nargin < 2
    dim = [];
    rephase = 0;
end 

    fid = fopen(filenm,'r');
    
    dat = fread(fid,inf,'int32');
    
    cplx = zeros([length(dat)/2 1]);
    
    
    for i=1:length(cplx)
       cplx(i) = dat(2*i-1) + sqrt(-1)*dat(2*i);
    end
    
    if isempty(dim);
        dim = length(cplx(:));
    end
    
    if prod(dim) ~= length(cplx(:))
        warning('dimensions do not match size of fid file...');
        return
    end
    
    %cplx = reshape(cplx,[dim(1)  dim(2)]);
    %cplx = permute(cplx, [1 3 2]);
   
    fclose(fid);
    
    if rephase
        % must rephase data
        
        % FT each slice
        for i=1:size(cplx, 3)
            
            cplxf = (fft2(cplx(:,:,i)));
            
            mag = abs(cplxf);
            phas = angle(cplxf);
            
           
            
            for j=2:2:size(cplxf,2); 
                for k=2:2:size(cplxf,1); 
                    
                    phas(k,j)=phas(k,j)+pi; 
                    
                    if(phas(k,j) > pi) 
                        phas(k,j) = phas(k,j) - 2*pi; 
                    end
                end 
            end;
            
            for j=1:2:size(cplxf,2); 
                for k=1:2:size(cplxf,1); 
                    
                    phas(k,j)=phas(k,j)+pi; 
                    
                    if(phas(k,j) > pi) 
                        phas(k,j) = phas(k,j) - 2*pi; 
                    end
                end 
            end;
            
            
            cplx(:,:,i) = (ifft2(mag .* exp(sqrt(-1) * phas)));
            
        end
        
                
    end
    