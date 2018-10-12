% Import needed slice or make a slice from imported data

for i=1:1:MaxFrames

    nameCur = strcat(filesbasement,num2str(Frame,'%04.f'),'.h5');

    attr = h5readatt(nameCur,'/Pid1','Parameters'); % Thread MASTER (0) always outputs its data
    gridno = attr(1);
    level = attr(2);
    mx = attr(3);
    my = attr(4);
    mz = attr(5);
    dx = attr(9);
    dy = attr(10);
    dz = attr(11);
    t = attr(12);
    iframe = attr(13);
    meqn = attr(14);
    lx = attr(15);
    ly = attr(16);
    lz = attr(17);

    fprintf('Working on data for Frame: %d at time: %d s for file: %s\n',Frame,t,nameCur);
   
    %-------------- Make meridional (y) slice ----------%

    if (strcmp(flagslice,'meridional'))

        if(isempty(listofproc))
        fprintf('Going through all cores to find threads with meridional slice of %d m ... \n',slicekm); 
        id = 1;
        for ii=1:1:lx
        for j=1:1:ly

        namedataset = strcat('/Pid',num2str(id));

        try
        a = hdf5read(nameCur,namedataset);

        attr = h5readatt(nameCur,namedataset,'Parameters');
        xlower = attr(6);
        ylower = attr(7);
        zlower = attr(8);
        xhigher = xlower + (mx-1)*dx;
        yhigher = ylower + (my-1)*dy;
        zhigher = zlower + (mz)*dz;

        if ((slicekm>=ylower) && (slicekm<=yhigher))
        namedataset;
        listofproc = [listofproc,id];
        myylower = ylower;
        myyhigher = yhigher;
        end

        catch
        fprintf('There is no data for thread:  %s \n',namedataset);
        end            
        id = id + 1;
        end
        end

        sliceinID = (slicekm-myylower)/dy + 1;

        if(isempty(listofproc) || isinf(sliceinID) || (floor(sliceinID) ~= sliceinID))
        error('There is no cell in y direction of %d m \n',slicekm);
        else
        %printf('Next cores contain meridional slice of %d m:\n',slicekm);
        listofproc
        end

        end

        datafullset = [];
        for ii=1:1:numel(listofproc)
        namedataset = strcat('/Pid',num2str(listofproc(ii)));
        namedataset;
        tempp = hdf5read(nameCur,namedataset);
        data = squeeze(tempp(:,sliceinID,:,:));
        datafullset = [datafullset;data];
        end
        pause
    end    

    %-------------- Make zonal (x) slice ----------%      

    if (strcmp(flagslice,'zonal'))

        if(isempty(listofproc))
        fprintf('Going through all cores to find threads with zonal slice of %d m ... \n',slicekm);

        id = 1;
        for ii=1:1:lx
        for j=1:1:ly

        namedataset = strcat('/Pid',num2str(id));

        try       
        a = hdf5read(nameCur,namedataset);
        attr = h5readatt(nameCur,namedataset,'Parameters');
        xlower = attr(6);
        ylower = attr(7);
        zlower = attr(8)  ; 
        xhigher = xlower + (mx-1)*dx;
        yhigher = ylower + (my-1)*dy;
        zhigher = zlower + (mz)*dz;

        if ((slicekm>=xlower) && (slicekm<xhigher))
        listofproc = [listofproc,id];
        myxlower = xlower;
        myxhigher = xhigher;
        end

        catch
        fprintf('There is no data for thread: %s \n',namedataset);
        end  
        id = id + 1;
        end
        end

        sliceinID = (slicekm-myxlower)/dx + 1;

        if(isempty(listofproc) || isinf(sliceinID) || (floor(sliceinID) ~= sliceinID))
        error('There is no cell in x direction of %d m \n',slicekm);
        else
        fprintf('Next cores contain zonal slice of %d m:\n',slicekm);
        listofproc
        end

        end

        datafullset = [];
        for ii=1:1:numel(listofproc)
        namedataset = strcat('/Pid',num2str(listofproc(ii)));
        namedataset;
        tempp = hdf5read(nameCur,namedataset);
        data = squeeze(tempp(sliceinID,:,:,:));
        datafullset = [datafullset;data];
        end

    end

    %-------------- Make horizontal slice --------------%      
    % Horizontal slicing takes most of time because     %
    % it requeres loading data from every processor     %
    %---------------------------------------------------%  

    if (strcmp(flagslice,'horizontal'))

        cellataltitude = mz*slicekm/(mz*dz);

        if(isinf(cellataltitude) || (floor(cellataltitude) ~= cellataltitude))
        error('Error: There is no cell at altitude of %d m \n',cellataltitude);
        end

        id = 1;
        dataset=[];
        datafullset=[];
        jone = 1;
        for ii=1:1:lx
        for jj=1:1:ly
        namedataset = strcat('/Pid',num2str(id));
        tempp = hdf5read(nameCur,namedataset);
        dataset = [dataset,squeeze(tempp(:,:,cellataltitude,:))];
        id = id+1;
        end

        datafullset = [datafullset;dataset];
        size(datafullset)
        dataset=[];
        end

    end

    %-------------- Plot horizontal slice --------------%      
    % This is to work with horizontal sliced output     %
    % it requeres loading data from every processor     %
    %---------------------------------------------------%  

    if (strcmp(flagslice,'horizontalslice'))

        cellataltitude = 1;

        id = 1;
        dataset=[];
        datafullset=[];
        jone = 1;
        for ii=1:1:lx
        for jj=1:1:ly
        namedataset = strcat('/Pid',num2str(id));
        tempp = hdf5read(nameCur,namedataset);
        dataset = [dataset,squeeze(tempp(:,:,cellataltitude,:))];
        id = id+1;
        end

        datafullset = [datafullset;dataset];
        size(datafullset);
        dataset=[];
        end
  
    end
    
    %-------------- Input full 3D domain  --------------%
    if (strcmp(flagslice,'full'))

        cellataltitude = 1;
        id = 1;
        dataset=[];
        datafullset=[];
        for ii=1:1:lx
        for jj=1:1:ly
        namedataset = strcat('/Pid',num2str(id));
        tempp = hdf5read(nameCur,namedataset);
        dataset = [dataset,squeeze(tempp(:,:,:,:))];
        id = id+1;
        end

        datafullset = [datafullset;dataset];
        size(datafullset);
        dataset=[];
        end
        
        fprintf('Full 3D domain is loaded\n');
        pause
    end
    
    
    
    if(figuresoutput)
    plotting
    if(videooutput)
    mov(it) = getframe(gcf);
    close all   
    it = it +1;
    end
    end

    Frame = Frame + 1;

end