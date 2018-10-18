
! =========================================================
      subroutine outslicever3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,&
     &     mz,xlower,ylower,zlower,dx,dy,dz,q,t,iframe)
! =========================================================

! Routine to save compressed chunked output from MAGIC !

! Please note that:
! The effectiveness of compression/writing/reading of hdf5 files
! depends on choice of chunk sizes. Some tradeoff should be done
! Current version uses a simple algorithm to define chunks sizes

     use mpi
     use hdf5
     implicit double precision (a-h,o-z)
    
     !------------------ HDF variables ------------------!
     integer(hid_t) :: plist_id      ! property list identifier
     integer(hid_t) :: dcpl          ! property list identifier
     integer(hid_t) :: file_id       ! file identifier
     integer(hid_t) :: dataset_id    ! dataset identifier
     integer(hid_t) :: dataspace_id  ! dataspace identifier
     double precision, allocatable :: attr_data(:,:)   ! data to write
     double precision, dimension(17) :: attr_datacur
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/17/) ! Attribute dimension
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
     INTEGER(HID_T) :: attr_id ! Attribute identifier
     INTEGER(HID_T) :: aspace_id ! Attribute dataspace identifier
     INTEGER(HID_T) :: atype_id ! Attribute dataspace identifier
     INTEGER :: arank = 1 ! Attribute rank
     !---------------------------------------------------! 

     !------------------ Filter variables ---------------!
     double precision, allocatable, target :: data(:,:,:,:)           ! data   
     integer(hsize_t), dimension(4) :: cdims = (/1,1,1,1/) ! chunks data dimensions
     ! integer, dimension(1:1) :: cd_values                  ! auxiliary data for the filter
     ! integer(size_t) :: nelmts                             ! number of elements in cd_values
     ! integer :: flags                                      ! bit vector specifying certain general properties of the filter
     ! integer(size_t) :: namelen = 180                      ! anticipated number of characters in name
     ! character(len=180) :: name                            ! name of the filter
     ! integer :: filter_id                                  ! filter identification number
     !---------------------------------------------------!

     !------------------ miscellaneous ------------------!
     character(len=3) :: c                                     ! dataset name for specific rank
     character(len=10) :: dataset_name
     integer :: rank = 4                                       ! data rank. q is 4D
     integer, allocatable :: idarray(:)
     character(mpi_max_processor_name) hostname
     dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,1-mbc:maxmz+mbc,meqn)
     integer(hsize_t), dimension(4) :: dimsf ! data dataset dimensions
     integer :: i,j,k,l,m,info,idd,arraysize,iddd
     character*20 fname
     common /mpicomm/ mpi_comm_3d, lx, ly, lz
     common /mpi_proc_info/ np, id
     !---------------------------------------------------!

     ! initialize HDF5 fortran interface
    call h5open_f(ierr)
    
    ! define size of q for every core
    dimsf(1) = mx
    dimsf(2) = my
    dimsf(3) = mz
    dimsf(4) = meqn
    
    if(id==0) then
    print*,dimsf(1)
    print*,dimsf(2)
    print*,dimsf(3)
    print*,dimsf(4)
    end if
    
    arraysize = 11
    
    allocate (data(dimsf(1),dimsf(2),dimsf(3),dimsf(4)))
    allocate (attr_data(3,arraysize))
    allocate (idarray(arraysize))
    
    idarray = (/0,12,36,60,84,108,132,156,180,204,228/)
      
     info = mpi_info_null
     fname = 'fort.qv' &
     & // char(ichar('0') + mod(iframe/1000,10)) &
     & // char(ichar('0') + mod(iframe/100,10)) &
     & // char(ichar('0') + mod(iframe/10,10)) &
     & // char(ichar('0') + mod(iframe,10)) &
     & // '.h5'
              
      ! Check data for very small values and 
      do k=1,mz
      do j=1,my
      do i=1,mx
      do m=1,meqn
      if (abs(q(i,j,k,m)) .lt. 1d-99) then
      q(i,j,k,m) = 0.d0
      end if
      data(i,j,k,m) = q(i,j,k,m)
      end do
      end do
      end do
      end do
          
     ! create datatype for the attribute
 	 ! Copy existing datatype
 	 ! H5T_NATIVE_CHARACTER - copy this datatype
     ! atype_id - copy datatype to this variable
	 call h5tcopy_f(h5t_native_double,atype_id,ierr)
    
    ! have id 0 creates hdf5 data layout and write all attributes
    if (id == 0) then
    
    ! This is a simple chunking of data
    if (dimsf(1) == 1) then
    cdims(1) = 1
    else if(mod(dimsf(1),5) == 0) then
    cdims(1) = dimsf(1)/5
    else if(mod(dimsf(1),2) == 0) then
    cdims(1) = dimsf(1)/2
    else
    cdims(1) = dimsf(1)
    end if
    
    if (dimsf(2) == 1) then
    cdims(2) = 1
    else if(mod(dimsf(2),5) == 0) then
    cdims(2) = dimsf(2)/5
    else if(mod(dimsf(2),2) == 0) then
    cdims(2) = dimsf(2)/2
    else
    cdims(2) = dimsf(2)
    end if
    
    if (dimsf(3) == 1) then
    cdims(3) = 1
    else if(mod(dimsf(3),5) == 0) then
    cdims(3) = dimsf(3)/5
    else if(mod(dimsf(3),2) == 0) then
    cdims(3) = dimsf(3)/2
    else
    cdims(3) = dimsf(3)
    end if
    
    if (meqn == 1) then
    cdims(4) = 1
    else if(mod(meqn,5) == 0) then
    cdims(4) = meqn/5
    else if(mod(meqn,2) == 0) then
    cdims(4) = meqn/2
    else
    cdims(4) = meqn
    end if
    
    ! Calculate xlower, ylower, zlower for every processor
     idd = 1
     iddd = 0
     k = 0 ! Current MAGIC does not support paralleling in vertical direction
     ! Only MASTER saves attributes, so it calculates xlow,ylow,zlow by itself
     do i=0,lx-1
     do j=0,ly-1
     !do kk=0,lz-1
     if (any(idarray==iddd)) then
        attr_data(1,idd) = i*mx*dx
        attr_data(2,idd) = j*my*dy
        attr_data(3,idd) = k*mz*dz
        idd = idd + 1
     end if
     iddd = iddd + 1
     !end do
     end do
     end do
     
     print*,'I AM DOING JOB1'
    
    ! create scalar dataspace for the attribute
    call h5screate_simple_f(arank,adims,aspace_id, ierr)
	
         ! create the hdf5 file
         ! filename: name of current file
         ! h5f_acc_trunc_f: rewrite if file exists
         ! file_id: identifier for file
         call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, ierr)
     
         ! create the dataspace for the dataset
         ! rank: rank of datasets
         ! dimsf: size of every dimension
         ! dataspace_id: identifier of dataspace
         call h5screate_simple_f(rank, dimsf, dataspace_id, ierr)
         
         ! create properties variable for the data
         ! h5p_dataset_create_f: property to create a data dataset
         ! dcpl: variable for data dataset properties
         call h5pcreate_f(h5p_dataset_create_f, dcpl, ierr)
        
         ! attribute the chunk size
         ! dcpl: link this property variable with chunk size
         ! 4: dimension of q
         ! cdims: dimensions of chunk in every direction
         call h5pset_chunk_f(dcpl, 4, cdims, ierr)
         
         ! attribute the compression type (GZIP compression)
         ! dcpl: link this property variable with filter
         ! 6: compression rank
!        call h5pset_deflate_f(dcpl, 6, ierr)

         ! attribute time of allocation of space for data in datasets
         ! h5d_alloc_time_early_f - allocate all space when the dataset is created
         call h5pset_alloc_time_f(dcpl, h5d_alloc_time_early_f, ierr)
         
        attr_datacur(1) = gridno
        attr_datacur(2) = level
        attr_datacur(3) = mx
        attr_datacur(4) = my
        attr_datacur(5) = mz
        attr_datacur(9) = dx
        attr_datacur(10) = dy
        attr_datacur(11) = dz
        attr_datacur(12) = t
        attr_datacur(13) = iframe
        attr_datacur(14) = meqn
        attr_datacur(15) = lx
        attr_datacur(16) = ly
        attr_datacur(17) = lz
        
        print*,'I AM DOING JOB1'
! create name for every dataset
   do i=1,arraysize

        attr_datacur(6) = attr_data(1,i) ! xlower
        attr_datacur(7) = attr_data(2,i) ! ylower
        attr_datacur(8) = attr_data(3,i) ! zlower

         write(c,"(i0)") idarray(i)+1
         dataset_name = "Pid" // trim(c)

         ! create dataset for this processor (based on id)
         ! file_id: name of file where to create dataset
         ! dataset_name: name of the dataset
         ! h5t_native_integer: type of data in dataset
         ! dataspace_id: dataspace used for dataset
         ! dataset_id: identifier for dataset
         ! dcpl_id: dataset creation property list
         call h5dcreate_f(file_id, dataset_name, h5t_native_double, &
                             dataspace_id, dataset_id, ierr, dcpl_id=dcpl)
                             
         ! Create attribute for dataset
         ! attr_id - ID of attribute
         call h5acreate_f(dataset_id,"Parameters",atype_id,aspace_id,attr_id, ierr)
         
         data_dims(1) = 17
         call h5awrite_f(attr_id, atype_id, attr_datacur, data_dims, ierr)
     
         ! close attribute
         call h5aclose_f(attr_id, ierr)                     

         ! close dataset
         call h5dclose_f(dataset_id, ierr)

   enddo
         print*,'I AM DOING JOB3'
         ! close access to the dataspace for attribute
         call h5sclose_f(aspace_id, ierr)
         
         ! close the dataspace
         call h5sclose_f(dataspace_id, ierr)
    
         ! close the properties variable
         call h5pclose_f(dcpl, ierr)
    
         ! close the file
         call h5fclose_f(file_id, ierr)
   end if
   
      ! mpi barrier to make sure everything is synched
      call mpi_barrier(mpi_comm_world, ierr)
     
      ! Now every processor is writing its own attributes and data to its dataset

      ! setup file access property variable with parallel i/o access
      ! plist_id: property variable
      ! comm - mpi communicator for mpi/io
      ! info - info regarding file access patterns and file system specifications
     call h5pcreate_f(h5p_file_access_f, plist_id, ierr)
     call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, info, ierr)

!      call h5pset_deflate_f(plist_id, 6, ierr)

     ! open hdf5 file for current time
     ! filename: filename of current hdf5 file
     ! h5f_acc_rdwr_f: open to read and write
     ! file_id: file identifier
     ! plist_id: file access property list
     call h5fopen_f(fname, h5f_acc_rdwr_f, file_id, ierr, plist_id)
  
     ! close the property list
     call h5pclose_f(plist_id, ierr)
   
     ! create the dataset names based on mpi rank
     write(c,"(i0)") id + 1
     dataset_name = "Pid" // trim(c)
     
     if (any(idarray==id)) then
     
     ! open dataset (each processor opens its own dataset)
     ! file_id: hdf5 file identifier
     ! dataset_name: dataset which belongs to this processor
     ! dataset_id: identifier for dataset
     call h5dopen_f(file_id, dataset_name, dataset_id, ierr)
        
     ! create properties variable
     ! h5p_dataset_xfer_f: property for raw data transfer
     call h5pcreate_f(h5p_dataset_xfer_f, plist_id, ierr)
     ! set collective mpio model
     ! h5fd_mpio_collective_f: collective is usually faster (OK to use it)
     ! call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, ierr)
     call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)
     
     ! write data to dataset
     ! dataset_id: identifier of dataset
     ! h5t_native_integer: type of data in memory which we want to write to file
     ! data: data by itself
     ! dimsf: dimensions of data we want to write to file
     ! xfer_prp = plist_id: data transfer property variable
     call h5dwrite_f(dataset_id, h5t_native_double, data, dimsf, ierr, xfer_prp = plist_id)
     
     print*,'we saved data!'
     print*,id
     
     call h5pclose_f(plist_id, ierr)
     call h5dclose_f(dataset_id,ierr)
      
    ! mpi barrier to make sure everything is synched
    ! call mpi_barrier(mpi_comm_world, ierr)
     end if

     call h5fclose_f(file_id, ierr)
     
      ! mpi barrier to make sure everything is synched
      ! call mpi_barrier(mpi_comm_world, ierr)
		
	  ! close fortran interface
      call h5close_f(ierr)
		
      return
      end