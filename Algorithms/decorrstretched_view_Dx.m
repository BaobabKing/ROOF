function decorrstretched_view_Dx(data)

   
    A = data;
    for k=1:size(data,3)
      A(:,:,k) = A(:,:,k) - mean2(A(:,:,k));
    end
    B = decorrstretch(A);  
    B = uint16(B);

   
    % For Axon Dx data from nCyte Prior platform 
    TOL=[0.001 0.999];
    for k=1:size(data,3)
        LH = stretchlim(B(:,:,k),TOL);
        LH(2) = max(LH(2),0.0155);  % normalized signal intensity must be greater than 1000/(2^16-1)=0.0153
        B(:,:,k) = imadjust(B(:,:,k),LH,[]);
    end
             
    RGB = cat(3,B(:,:,3)+B(:,:,2), B(:,:,4)+B(:,:,2), B(:,:,1)+B(:,:,2));
    figure, imshow(RGB); truesize; axis on; impixelinfo;
       
   
  
