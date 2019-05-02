 Module msd
     contains
     Subroutine calc_msd(xarray,yarray,zarray,Natoms,Nframes)
         Implicit none
         real(kind=8),dimension(:),intent(in)::xarray,yarray,zarray
         integer,intent(in)::Natoms,Nframes
         real(kind=8),dimension(3)::xyz
         real(kind=8),dimension(:,:,:),allocatable::trajectory
         integer::i,j,k,Nall,Ni, nf,na, iframe,jframe
         real(kind=8),dimension(:),allocatable::msd
         real(kind=8)::dt
         character(LEN=80)::outfilename
         Nall=Natoms*(Nframes)
         nf=1
         na=1
         write(*,*), 'nframes', Nframes
         ALLOCATE(trajectory(3,Natoms,Nframes+1))
         Do Ni =1,Nall
             xyz(1)=(xarray(Ni))
             xyz(2)=(yarray(Ni))
             xyz(3)=(zarray(Ni))
             trajectory(:,na,nf)=xyz(:)
             na = na + 1
             if (na .gt. Natoms) then
                 na = 1   
                 nf = nf + 1
             end if 
         end do
         ALLOCATE(msd(0:Nframes-1))
         write(*,*) "Computing mean squared displacement ( * = 10% )"
         msd(:)=0
         DO iframe=1,Nframes
             DO jframe=iframe,Nframes
                 DO i=1,Natoms
                     msd(jframe-iframe)=msd(jframe-iframe)+&
                     &SUM((trajectory(:,i,iframe)-trajectory(:,i,jframe))**2)
                 ENDDO
             ENDDO
         ENDDO
         ! normalize and write
         outfilename="MSDJK.dat"
         dt=0.5
         write(*,*) "Writing msd to ", outfilename
         open(11,FILE=outfilename)
         write(11,*) "# time [units as in input, e.g. fs],  MSD in Angstrom**2"
         DO iframe=0,Nframes-1
             jframe=Nframes-iframe
             msd(iframe)=(msd(iframe)/Natoms)/jframe
             write(11,*) iframe*dt,msd(iframe)
         ENDDO
         DEALLOCATE(msd)
 
 !        write(*,*) 'xyzarray', size(xarray), Natoms, Nframes, Nall
 !        write(*,*), 'trajectory', size(trajectory,1), size(trajectory,2), size(trajectory,3)
     End Subroutine
  End Module
