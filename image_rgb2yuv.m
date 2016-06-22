function str = image_rgb2yuv(im)

im=double(im);


    ar = im(:,:,1);
    ag = im(:,:,2);
    ab = im(:,:,3);


        ay=ar.*0.3+ag.*0.59+ab.*0.11;
        au=ar.*(-0.15)+ag.*(-0.29)+ab.*0.44;
        av=ar.*0.61+ag.*(-0.52)+ab.*(-0.096);


     im1 = ay;
     im2 = au;%(1:2:end,1:2:end)+128;
     im3 = av;%(1:2:end,1:2:end)+128;
   
%  im1=uint8(im1);
%  im2=uint8(im2);
%  im3=uint8(im3);
str(:,:,1)=im1;
str(:,:,2)=im2;
str(:,:,3)=im3;


end

