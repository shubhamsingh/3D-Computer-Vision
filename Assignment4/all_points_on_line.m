%Reference mathworks.com
function [x_coord,y_coord,marker] = all_points_on_line(point)
if (abs(point(4)-point(2)) > abs(point(3)-point(1)))                                      
    x0 = point(2);y0 = point(1); x1 = point(4);y1=point(3); 
    marker =1;                                              
else
    x0 = point(1);y0 = point(2); x1 = point(3);y1=point(4);
    marker = 0; 
end
if(x0 >x1)
    temp1 = x0; x0 = x1; x1 = temp1;
    temp2 = y0; y0 = y1; y1 = temp2;
end
dx = abs(x1 - x0) ;                              
dy = abs(y1 - y0);                               
sx = sign(x1 - x0);                              
sy = sign(y1 - y0);                              

x = x0; y = y0;                                  
param = 2*dy - dx ;                              
for i = 0:dx-1                                   
    x_coord(i+1) = x;                            
    y_coord(i+1) = y;
    param = param + 2*dy;                        
    if (param >0)                                
        y = y +1*sy;                             
        param = param - 2*(dx );                        
    end
    x = x + 1*sx;                                
end
end