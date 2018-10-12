function calcsliceh5(slice,dx,dy,dz,lx,ly,lz,mx,my,mz,flagslice)

id = 0; % do not change this
f = 0; % do not change this

fprintf('Output these threads: \n')

if (flagslice=='y')
for ii=0:1:lx-1
for jj=0:1:ly-1
xlower = ii*mx*dx;
ylower = jj*my*dy;
xhigher = xlower + (mx)*dx;
yhigher = ylower + (my)*dy;
if ((slice>=ylower) && (slice<yhigher))
f = f+1;
id
end
id = id+1;
end
end
fprintf('Total number of cores containing this slice is: %d \n',f)
end

id = 0; % do not change this
f = 0; % do not change this

if (flagslice=='x')
for ii=0:1:lx-1
for jj=0:1:ly-1
xlower = ii*mx*dx;
ylower = jj*my*dy;
xhigher = xlower + (mx)*dx;
yhigher = ylower + (my)*dy;
if ((slice>=xlower) && (slice<xhigher))
f = f+1;
id
end
id = id+1;
end
end
fprintf('Total number of cores containing this slice is: %d \n',f)
end

end

