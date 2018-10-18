c
c
c =========================================================
      subroutine out3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &                 xlower,ylower,zlower,dx,dy,dz,q,t,iframe)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 3 dimensions to a hdf file as a scientific data set.
c     # See http://hdf.ncsa.uiuc.edu/ for more info on HDF.
c     # The results in this output file can be plotted in MATLAB
c     # using the "plotclaw3" script.
c
c     # Revised 2003 by Peter Blossey
c     # Adapted code written by Sorin Mitran to standard F77 CLAWPACK style
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'
c
      parameter   (nDim = 3)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      character*14 fname
      character*13 qname2
      character*9  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
      integer    sd_id, sds_id, sds_rank, sds_dims, sds_start,
     $     sds_stride, sds_edges
      dimension  mstart(nDim), mtotal(nDim), sds_dims(nDim),
     $     sds_stride(nDim), sds_edges(nDim), sds_start(nDim)

c     # HDF: x- and y- dimensions are reversed when output to HDF file.
      dimension  qbuf(21), qsend(mz,my)
c
c     # HDF: Declare external HDF functions
c
      integer  sfstart, sfcreate, sfwdata, sfendacc, sfend
      external sfstart, sfcreate, sfwdata, sfendacc, sfend
c
c     # HDF: Set up HDF constants
      integer 	DFACC_READ, DFACC_WRITE, DFACC_CREATE 
      parameter(DFACC_READ = 1, DFACC_WRITE = 2, DFACC_CREATE = 4)

      integer   DFNT_FLOAT64, DFNT_INT32
      parameter(DFNT_FLOAT64 = 6, DFNT_INT32 = 24)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
      common /mpicomm/ mpi_comm_3d, lx, ly, lz, mtotal, mstart
      common /mpi_proc_info/ np, id
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # MPI: Let processor zero create the HDF file and write the grid
C     parameters to the HDF file.  
c
      if (id.eq.0) then
c     # first create the file name and open file
c     
         fname = 'fort.q'
     &        // char(ichar('0') + mod(iframe/1000,10)) 
     &        // char(ichar('0') + mod(iframe/100,10)) 
     &        // char(ichar('0') + mod(iframe/10,10)) 
     &        // char(ichar('0') + mod(iframe,10))
     &        // '.hdf'
c     
c     # Specify grid number and create a string which will describe this
c     # grid in HDF file.  This could be an input for simulations with
c     # multiple grids, as in AMRCLAW.  
c     
         ngrids_out=1
         qname = 'grid_'
     &        // char(ichar('0') + mod(ngrids_out/1000,10)) 
     &        // char(ichar('0') + mod(ngrids_out/100,10)) 
     &        // char(ichar('0') + mod(ngrids_out/10,10)) 
     &        // char(ichar('0') + mod(ngrids_out,10))
c     
c     # HDF: create hdf file.
         sd_id = sfstart(fname, DFACC_CREATE)
         if (sd_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to create ', fname,' (call to sfstart)'
            STOP
         end if
c     
c     # HDF: create a data set for parameters describing q in HDF file.
c
         sds_rank = 1
         sds_dims(1) = 21
      
         sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,sds_rank,sds_dims)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 
     &           'Failed to create scientific data set in HDF file'
            STOP
         end if
c
c     # HDF: set up parameters describing data set.
c
         sds_start(1)  = 0
         sds_edges  = sds_dims(1)
         sds_stride = 1        

         qbuf(1) = ngrids_out
         qbuf(2) = nDim
         qbuf(3) = t
         qbuf(4) = meqn
         qbuf(5) = 1.
         qbuf(6) = mtotal(1)
         qbuf(7) = mtotal(2)
         qbuf(8) = mtotal(3)
         qbuf(9) = 0.
         qbuf(10) = xlower
         qbuf(11) = ylower
         qbuf(12) = zlower
         qbuf(13) = 0.
         qbuf(14) = xlower+mtotal(1)*dx
         qbuf(15) = ylower+mtotal(2)*dy
         qbuf(16) = zlower+mtotal(3)*dz
         qbuf(17) = 0.
         qbuf(18) = dx
         qbuf(19) = dy
         qbuf(20) = dz
         qbuf(21) = 0.
         istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qbuf)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to write grid data (call to sfwdata)'
            STOP
         end if
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close grid data (call to sfendacc)'
            STOP
         end if
c
      end if
c
c     # MPI: Wait for grid information to be written.
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c     
c     # Loop over fields in q
c     
      do m = 1,meqn
         
         if (id.eq.0) then
c     
c     # HDF: create a data set for parameters describing q in HDF file.
            qname2 = qname // '_'
     &           // char(ichar('0') + mod(m/100,10)) 
     &           // char(ichar('0') + mod(m/10,10)) 
     &           // char(ichar('0') + mod(m,10))
c
c     # HDF: Reverse dimensions when storing arrays because of different
c     #      conventions between c and FORTRAN as to which dimension should
c     #      be stored first.  Reversing the dimensions here will make
c     #      the x-direction first when reading into MATLAB.
c
            sds_rank = nDim
            sds_dims(1) = mtotal(3)
            sds_dims(2) = mtotal(2)
            sds_dims(3) = mtotal(1)
c     
c     # HDF: Create data set for this component of q.
            sds_id = sfcreate(sd_id, qname, DFNT_FLOAT64, sds_rank,
     &           sds_dims)
            if (sds_id.eq.FAIL) THEN
               WRITE(*,*) 'Failed to create data set number ', m,
     &              ' in HDF file'
               STOP
            end if
c
         end if
c
         do i = 1,mtotal(1)
c
c     # If this processor has data at x-slice i, load into qsend
            if (((mstart(1)+1).le.i).and.((mstart(1)+mx).ge.i)) then
               do j = 1,my
                  do k = 1,mz
                     qsend(k,j) = q(i-mstart(1),j,k,m)
                  end do
               end do
            end if

            if (id.eq.0) then
c     # If processor 0, call routine to receive and write this slice of data.
               call gather_slice(lz,ly,mtotal(3),mtotal(2),mz,my,i,
     &              sds_id,qsend)
            elseif (((mstart(1)+1).le.i).and.((mstart(1)+mx).ge.i)) then
c     # If this processor has data at x-slice i, then send to processor 0.
               call mpi_send(qsend, my*mz, MPI_DOUBLE_PRECISION,
     &              0, id, MPI_COMM_WORLD, ierr)
            end if
c
c     # MPI: Wait for this slice is written to disk.
            call mpi_barrier(MPI_COMM_WORLD,ierr)
c
         end do
c
c     # HDF: Close the data set
         if (id.eq.0) then
            istat = sfendacc(sds_id)  
            if (istat.eq.FAIL) then
               WRITE(*,*) 'Failed to close SDS number ', m,
     &              ' (call to sfendacc)'
               STOP
            end if
         end if
c
      end do
c     
c     # HDF: Close HDF file.
      if (id.eq.0) then
         istat = sfend(sd_id)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close ',fname,' (call to sfend)'
            STOP
         end if
      end if

      return
      end 

c     #################################################################
      subroutine gather_slice(nz,ny,mztot,mytot,mz,my,i,sds_id,qsend)

      implicit none
      include 'mpif.h'
c
c     # Note that nz and ny are the dimensions of the processor array
c     #           mz and my are the dimensions of the array on this processor
c     #           mztot and mytot are the overall dimensions.
      integer  i, my, mz, mytot, mztot, nz, ny, sds_id
      integer  mprocid(nz,ny), ireq(nz,ny)
      double precision  qout(mztot,mytot), qsend(mz,my)
      double precision  qgather((mz+1)*(my+1),nz,ny)
c
      integer    idcomm, ii, jj, kk, j, k, icoord(3)
      integer    index, ierr
      integer    joffset, koffset, jlength, klength, nrecv, nDim
      parameter  (nDim = 3)
c
c     # HDF: Declare variables that describe datasets and HDF files.
      integer    sds_start(nDim), sds_edges(nDim), sds_stride(nDim)
      integer   istat, istatus(MPI_STATUS_SIZE)
c
c     # HDF: Declare external HDF functions
      integer  sfwdata
      external sfwdata
c
c     # HDF: Set up HDF constants
      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c     # MPI: MPI-Related common blocks.
      integer  mpi_comm_3d, lx, ly, lz
      integer  mstart(nDim), mtotal(nDim)
      common /mpicomm/ mpi_comm_3d, lx, ly, lz, mtotal, mstart

      integer  maxprocs
      parameter (maxprocs = 360)
      integer  mproc_start(0:maxprocs-1,nDim), mx_start(maxprocs)
      integer  mproc_length(0:maxprocs-1,nDim), mx_length(maxprocs)
      common /mpi_length_comm/ mproc_start,mproc_length,
     &     mx_start,mx_length

      integer  np, id
      common /mpi_proc_info/ np, id
c
c     # MPI: Post receive from various processors that contain a piece
C     # of x-slice i.
      do ii = 1,lx
         if ((i.ge.(mx_start(ii)+1)).and.
     &        (i.le.(mx_start(ii)+mx_length(ii)))) icoord(1) = ii-1
      end do
      do j = 1,ly
         do k = 1,lz
            icoord(2) = j-1
            icoord(3) = k-1
            call mpi_cart_rank(mpi_comm_3d,icoord,idcomm,ierr)
            nrecv = mproc_length(idcomm,2)*mproc_length(idcomm,3)
            call mpi_irecv(qgather(1,k,j), nrecv, MPI_DOUBLE_PRECISION,
     &           idcomm, idcomm, MPI_COMM_WORLD, ireq(k,j), ierr)
            mprocid(k,j) = idcomm
         end do
      end do
c
c     # If this processor has data at x-slice i, send to itself.
      if (((mstart(1)+1).le.i).and.
     &     ((mstart(1)+mproc_length(id,1)).ge.i)) then
         call mpi_send(qsend, my*mz, MPI_DOUBLE_PRECISION,
     &        0, id, MPI_COMM_WORLD, ierr)
      end if
c
c     # MPI: Wait for message completion.  Then, sort contribution to
c     # qgather from various processors into qout for output to HDF file.
      do j = 1,ly
         do k = 1,lz
            call mpi_wait(ireq(k,j),istatus,ierr)
            idcomm = mprocid(k,j)
            joffset = mproc_start(idcomm,2)
            koffset = mproc_start(idcomm,3)
            jlength = mproc_length(idcomm,2)
            Klength = mproc_length(idcomm,3)
            index = 1
            do jj = 1,jlength
               do kk = 1,klength
                  qout(koffset+kk,joffset+jj) = qgather(index,k,j)
                  index = index + 1
               end do
            end do
         end do
      end do
c     
c     # HDF: set up parameters describing data set.
      sds_start(1)  = 0
      sds_start(2)  = 0
      sds_start(3)  = i-1

      sds_edges(1)  = mtotal(3)
      sds_edges(2)  = mtotal(2)
      sds_edges(3)  = 1

      sds_stride(1) = 1        
      sds_stride(2) = 1        
      sds_stride(3) = 1        
c     
c     # HDF: write current processor's slab of data to hdf file.
      istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qout)
      if (istat.eq.FAIL) then
         WRITE(*,*) 'Failed to write SDS (call to sfwdata)'
         STOP
      end if

      return
      end
