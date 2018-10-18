c     
c     =====================================================
      subroutine restart(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &     xlower,ylower,zlower,dx,dy,dz,q)
c     =====================================================
c     
c     # Initialize q using values from an old HDF output file.
c     # Copy the HDF output file to restart.data.hdf
c     # and call this routine from qinit.
c     
c     # See http://hdf.ncsa.uiuc.edu/ for more info on HDF.
c     
c     # Written 2003 by Peter Blossey 
c     
      implicit double precision (a-h,o-z)
      include 'mpif.h'
c     
      parameter   (nDim = 3)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &     1-mbc:maxmz+mbc, meqn)
      character*14 fname
c     
c     # HDF: Declare variables that describe datasets and HDF files.
c     
      integer    sd_id, sds_id, sds_start, sds_edges, sds_stride
      dimension  sds_start(nDim), sds_edges(nDim), sds_stride(nDim) 

c     # HDF: x- and y- dimensions are reversed when output to HDF file.
      dimension  qbuf(21), qout(mz,my,mx)
c     
c     # HDF: Declare external HDF functions
c     
      integer  sfstart, sfselect, sfrdata, sfendacc, sfend
      external sfstart, sfselect, sfrdata, sfendacc, sfend
c     
c     # HDF: Set up HDF constants
c     
      integer 	DFACC_READ
      parameter(DFACC_READ = 1)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
      logical   outt0
      common /restrt_block/ tinitial, iframe, outt0
c
c     # MPI: MPI params
c
      dimension  id_coords(3), mstart(3), mtotal(3), idimlength(3)
      dimension  isds_id(meqn), istatus(MPI_STATUS_SIZE)
      common /mpicomm/ mpi_comm_3d, lx, ly, lz, mtotal, mstart
      common /mpi_proc_info/ np, id
c     
c     # first create the file name and open file
c     
      fname = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
     &     // '.hdf'
      write(*,*) 'Restarting from ', fname
c
c     # MPI: get info about this processor's position in array of processors.
c
      call mpi_cart_coords(mpi_comm_3d,id,3,id_coords,ierr)
c
c     # HDF: Set up index, which selects field from HDF output file.
c     #      Used in calls to sfselect.  (First field in file has index=0.)
c
      index = 0
c     
c     # HDF: open hdf restart file.
c     
      if (mstart(1)+mstart(2)+mstart(3).eq.0) then
c
         sd_id = sfstart(fname,DFACC_READ)
         if (sd_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to open HDF file',
     &           ' (call to sfstart in restart3_mpi_hdf.f)'
            STOP
         end if
c     
c     # HDF: Select grid parameter dataset in HDF file.
c     
         sds_id = sfselect(sd_id,index)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to select data set for  variable ', 
     &           index, ' in restart HDF file',
     &                 ' on processor ', id
            WRITE(*,*) '(call to sfselect in restrt3_mpi_hdf.f)'
            STOP
         end if
c     
c     # HDF: Set up dimensions for double vector.
c     
         sds_start(1) = 0
         sds_edges(1) = 21
         sds_stride(1) = 1
c     
c     # HDF: read double vector from hdf file.
c     
         istat = sfrdata(sds_id,sds_start,sds_stride,sds_edges,qbuf)  
         if (istat.eq.FAIL) THEN
            WRITE(*,*) 'Failed to read variable ', index,
     &           ' from restart HDF file'
            WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
            STOP
         end if
c     
c     # HDF: End access to double vector in HDF file.
c     
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) THEN
            WRITE(*,*) 'Failed to end access to variable ', index,
     &           ' in restart HDF file'
            WRITE(*,*) '(call to sfendacc in restrt_hdf.f)'
            STOP
         end if
c
c     # HDF: Close HDF file.
c     
         istat = sfend(sd_id)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close ',fname,' (call to sfend)'
            STOP
         end if
c     
c     # Check parameters against those read in from claw3ez.data
c     
         mx_in = qbuf(6)
         my_in = qbuf(7)
         mz_in = qbuf(8)
         if ((mx_in .ne. mtotal(1)) .or. 
     &        (my_in .ne. mtotal(2)) .or. 
     &        (mz_in .ne. mtotal(3))) then
            stop 'rstart.f : grid dimensions not compatible'
         endif
c     
         meqn_in = qbuf(4)
         if (meqn_in .ne. meqn) then
            stop 'rstart.f : meqn not compatible'
         endif
c
c     # Read in new starting time.
c
         tinitial = qbuf(3)
c
         write(*,*) 'Restarting from file ', fname, ' at time ', 
     &        tinitial
         write(*,*) 'Initial condition will not be written to a ',
     &        'MATLAB output file'
         write(*,*)
c
      end if
c
c     # MPI: Wait for master process to finish reading in grid data.
c
      call mpi_barrier(mpi_comm_world,ierr)
c
c     # Since we are starting from an old field, don't write out initial
c     # conditions since it will just be a copy of the old file.
c
      outt0 = .false.
c     
c     # Loop over fields in q
c     
      do m = 1,meqn
c
c     # HDF : increment index.
c
         index = index + 1
c
c     # MPI : Loop over all processors to read in data.
c
         do n = 0,np-1
c
            if (id.eq.n) then
c     
c     # HDF: Open HDF file.
c     
               sd_id = sfstart(fname,DFACC_READ)
               if (sd_id.eq.FAIL) THEN
                  WRITE(*,*) 'Failed to open HDF file',
     &                 ' (call to sfstart in restart3_mpi_hdf.f)'
                  STOP
               end if
c     
c     # HDF: Select dataset in HDF file.
c     
               sds_id = sfselect(sd_id,index)
               if (sds_id.eq.FAIL) THEN
                  WRITE(*,*) 'Failed to select data set number ', 
     &                 index, ' in restart HDF file'
                  WRITE(*,*) '(call to sfselect in restrt3_mpi_hdf.f)',
     &                 ' on processor ', id
                  STOP
               end if
c     
c     # HDF: Set up dimensions for double array.
c     
               sds_start(1) = mstart(3)
               sds_start(2) = mstart(2)
               sds_start(3) = mstart(1)

               sds_edges(1) = mz
               sds_edges(2) = my
               sds_edges(3) = mx

               sds_stride(1) = 1
               sds_stride(2) = 1
               sds_stride(3) = 1
c     
c     # HDF: read double array from hdf file.
c     
               istat = sfrdata(sds_id,sds_start,sds_stride,sds_edges,
     &              qout)  
               if (istat.eq.FAIL) THEN
                  WRITE(*,*) 'Failed to read variable ', index,
     &                 ' from restart HDF file'
                  WRITE(*,*) '(call to sfrdata in restrt3_mpi_hdf.f)',
     &                 ' on processor ', id
                  STOP
               end if
c     
c     # HDF: End access to double array in HDF file.
c     
               istat = sfendacc(sds_id)  
               if (istat.eq.FAIL) THEN
                  WRITE(*,*) 'Failed to end access to variable ', index,
     &                 ' in restart HDF file'
                  WRITE(*,*) '(call to sfendacc in restrt_hdf.f)'
                  STOP
               end if
c
c     # HDF: Close HDF file.
c     
               istat = sfend(sd_id)
               if (istat.eq.FAIL) then
                  WRITE(*,*) 'Failed to close ',fname,' (call to sfend)'
                  STOP
               end if
c     
c     # Put the data into q(:,:,:,m).
c
               do k = 1,mz
                  do j = 1,my
                     do i = 1,mx
                        q(i,j,k,m) = qout(k,j,i)
                     end do
                  end do
               end do
c
            end if
c
c     # MPI: Sync up processors after each one reads its chunk of data
c
            call mpi_bcast(n,1,MPI_INTEGER,n,MPI_COMM_WORLD,ierr)
c
c     # MPI: end loop over processors.
c
         end do
c
c     # end loop over variables.
c
      end do

      return
      end
