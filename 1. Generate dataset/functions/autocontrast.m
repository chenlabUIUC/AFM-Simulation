% AUTOCONTRAST  Automatically adjusts contrast of images to optimum level.
function I = autocontrast(I0, l1, l2)
if nargin<2
    l1 = 0.02;
    l2 = 0.98;
elseif nargin<3
    l2 = 1-l1;
end
low_limit=l1;
up_limit=l2;
img=I0;
[m1 n1 r1]=size(img);
% img=double(img);
%--------------------calculation of vmin and vmax----------------------
for k=1:r1
    arr=sort(reshape(img(:,:,k),m1*n1,1));
    if low_limit == 0
        v_min(k) = 1;
    else
        v_min(k)=arr(ceil(low_limit*m1*n1));
    end
    v_max(k)=arr(ceil(up_limit*m1*n1));
end
%----------------------------------------------------------------------
if r1==3
    v_min=rgb2ntsc(v_min);
    v_max=rgb2ntsc(v_max);
end
%----------------------------------------------------------------------
img=(img-v_min(1))/(v_max(1)-v_min(1));
img(img>1) = 1;
img(img<0) = 0;
I = img;
end