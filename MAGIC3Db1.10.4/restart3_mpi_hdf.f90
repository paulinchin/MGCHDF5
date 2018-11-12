      subroutine restart(maxmx,maxmy,maxmz,meqn,mbc,mx,&
     &     my,mz,xlower,ylower,zlower,dx,dy,dz,q)
     
      use mpi
      use hdf5   
      implicit double precision (a-h,o-z)

     !------------------ HDF variables ------------------!
     integer(hid_t) :: plist_id      ! property list identifier
     integer(hid_t) :: dcpl          ! property list identifier
     integer(hid_t) :: file_id       ! file identifier
     integer(hid_t) :: dataset_id    ! dataset identifier
     integer(hid_t) :: dataspace_id  ! dataspace identifier
     double precision, dimension(17) :: attr_dataread
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/17/) ! Attribute dimension
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
     INTEGER(HID_T) :: attr_id ! Attribute identifier
     INTEGER(HID_T) :: aspace_id ! Attribute dataspace identifier
     INTEGER(HID_T) :: atype_id ! Attribute dataspace identifier
     INTEGER :: arank = 1 ! Attribute rank
     !---------------------------------------------------! 

     !------------------ Filter variables ---------------!
     double precision, allocatable, target :: data_out(:,:,:,:)           ! data   
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
 !    character(mpi_max_processor_name) hostname
 !    dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,1-mbc:maxmz+mbc,meqn)
     dimension q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc)
     integer(hsize_t), dimension(4) :: dimsf ! data dataset dimensions
     integer :: i,j,k,l,m,info,idd
     character*20 fname
     common /mpicomm/ mpi_comm_3d, lx, ly, lz
     common /mpi_proc_info/ np, id
     logical outt0
     common /restrt_block/ tinitial, iframe, outt0
     !---------------------------------------------------!

!     MPI: get info about this processor position in array of processors.
      call mpi_cart_coords(mpi_comm_3d,id,3,id_coords,ierr)     
      
     ! create the file name and open file
     ! initialize HDF5 fortran interface
     
    call h5open_f(ierr)
     
    info = mpi_info_null
    
    ! define size of q for every core
    dimsf(1) = meqn
    dimsf(2) = mx
    dimsf(3) = my
    dimsf(4) = mz     
    
    allocate (data_out(dimsf(1),dimsf(2),dimsf(3),dimsf(4)))
    
!	Set name for the file you start from

    fname = 'fort.q0600.h5'
!    'fort.q' &
!     &     // char(ichar('0') + mod(iframe/1000,10)) &
!     &     // char(ichar('0') + mod(iframe/100,10)) &
!     &     // char(ichar('0') + mod(iframe/10,10)) &
!     &     // char(ichar('0') + mod(iframe,10)) &
!     &     // '.h5'

!    write(*,*) 'Restart from ', fname
	
      ! setup file access property variable with parallel i/o access
      ! plist_id: property variable
      ! comm - mpi communicator for mpi/io
      ! info - info regarding file access patterns and file system specifications
     call h5pcreate_f(h5p_file_access_f, plist_id, ierr)
     call h5pset_fapl_mpio_f(plist_id, mpi_comm_3d, info, ierr)

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
     
     ! print*,dataset_name
     
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
     call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, ierr)
     
     ! write data to dataset
     ! dataset_id: identifier of dataset
     ! h5t_native_integer: type of data in memory which we want to write to file
     ! data_out: data by itself
     ! dimsf: dimensions of data we want to write to file
     ! xfer_prp = plist_id: data transfer property variable
     call h5dread_f(dataset_id, h5t_native_double, data_out, dimsf, ierr, xfer_prp = plist_id)

	 ! Put the data into q by every processor
     do m = 1,meqn
     do k = 1,mz
     do j = 1,my
     do i = 1,mx
     q(m,i,j,k) = data_out(m,i,j,k)
     end do
     end do
     end do
     end do
           
     call h5pclose_f(plist_id, ierr)
     call h5dclose_f(dataset_id,ierr)
     call h5fclose_f(file_id, ierr)
     
     deallocate(data_out)
  
     ! close fortran interface
     call h5close_f(ierr)
     
      return
      end
