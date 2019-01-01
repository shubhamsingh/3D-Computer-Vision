function [ou1] = warpping(in1,H)
   % Interpolation='linear'; Boundary='zero'
   % Reference https://in.mathworks.com/
   [x,y]=ndgrid(0:size(in1,1)-1,0:size(in1,2)-1);
   xd=x;
   yd=y;
   Tlocalx = mean(1) + H.T(1,1) * xd + H.T(1,2) *yd + H.T(1,3) * 1;
   Tlocaly = mean(2) + H.T(2,1) * xd + H.T(2,2) *yd + H.T(2,3) * 1;
   
   in1=double(in1);
   ImageSize=[size(in1,1) size(in1,2)];
   if(ndims(in1)==2), lo=1; else lo=3; end

    xBas0=floor(Tlocalx);
    yBas0=floor(Tlocaly);
    xBas1=xBas0+1;
    yBas1=yBas0+1;

    tx=Tlocalx-xBas0;
    ty=Tlocaly-yBas0;
    perc0=(1-tx).*(1-ty);
    perc1=(1-tx).*ty;
    perc2=tx.*(1-ty);
    perc3=tx.*ty;

    check_xBas0=(xBas0<0)|(xBas0>(size(in1,1)-1));
    check_yBas0=(yBas0<0)|(yBas0>(size(in1,2)-1));
    check_xBas1=(xBas1<0)|(xBas1>(size(in1,1)-1));
    check_yBas1=(yBas1<0)|(yBas1>(size(in1,2)-1));
    xBas0=min(max(xBas0,0),size(in1,1)-1);
    yBas0=min(max(yBas0,0),size(in1,2)-1);
    xBas1=min(max(xBas1,0),size(in1,1)-1);
    yBas1=min(max(yBas1,0),size(in1,2)-1);

    ou1=zeros([ImageSize(1:2) lo]);
    for i=1:lo
        Iin_one=in1(:,:,i);
        intensity_xyz0=Iin_one(1+xBas0+yBas0*size(in1,1));
        intensity_xyz1=Iin_one(1+xBas0+yBas1*size(in1,1));
        intensity_xyz2=Iin_one(1+xBas1+yBas0*size(in1,1));
        intensity_xyz3=Iin_one(1+xBas1+yBas1*size(in1,1));     
        intensity_xyz0(check_xBas0|check_yBas0)=0;
        intensity_xyz1(check_xBas0|check_yBas1)=0;
        intensity_xyz2(check_xBas1|check_yBas0)=0;
        intensity_xyz3(check_xBas1|check_yBas1)=0;
        Iout_one=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3;
        ou1(:,:,i)=reshape(Iout_one, ImageSize);
    end
end