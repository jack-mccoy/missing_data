! homemade EM to estimate mean and variance subject to missing data
! Fortran version

subroutine mvn_emf(Ey, estE, estR, NAindex, tol, maxiter, K, N)
	
implicit none

!declare input arguments	
integer K,N
double precision Ey(N,K), estE(K), estR(K,K)
logical NAindex(N,K)
double precision tol
integer maxiter

!declare other stuff
integer i,j,iter,NAcount,info
double precision Enew(K), Rnew(K,K)
double precision Eyy(K,K), thisEy(K), thisEyy(K,K)
double precision Rdist, Edist
logical thisNA(K), lastNA(K), thisNA2(K,K), thisDataNA(K,K), thisData2(K,K)
double precision, allocatable :: E1(:), E2(:)
double precision, allocatable :: R11(:,:), R12(:,:), R22(:,:), invR11R12(:,:)

!calculations start here
write(*,'(A)') 'Starting EM:'

!allocate first, because we always deallocate first in the loop
allocate(E1(1),E2(1),R11(1,1),R12(1,1),R22(1,1),invR11R12(1,1)) 

iterloop: do iter=1,maxiter

   Eyy=0
   lastNA=.false.
   
   
   !compute means and variances for every observation (E-step)
   iloop: do i=1,N

      NAcount=count(NAindex(i,:))
      if (NAcount == K) then !no data

         thisEy=estE
         thisEyy=estR

      elseif (NAcount == 0) then !complete data

         thisEy=Ey(i,:)
         thisEyy=0

      else !some data

         thisNA=NAindex(i,:)
         if (any(thisNA .neqv. lastNA)) then !can save the below computations if missingness pattern is the same
            !set up indices            
            do j=1,K
               thisNA2(j,:)=thisNA(j) .and. thisNA
               thisDataNA(j,:)= (.not.thisNA(j)) .and. thisNA
               thisData2(j,:)= (.not.thisNA(j)) .and. (.not.thisNA)
            end do

            !re-allocatae
            deallocate(E1,E2,R11,R12,R22,invR11R12)
            allocate(E1(K-NAcount),E2(NAcount))
            allocate(R11(K-NAcount,K-NAcount))
            allocate(R12(K-NAcount,NAcount))
            allocate(R22(NAcount,NAcount))
            allocate(invR11R12(K-NAcount,NAcount))
            
            !calculate the key matrices
            R11=reshape(pack(estR,thisData2),(/K-NAcount,K-NAcount/))
            R12=reshape(pack(estR,thisDataNA),(/K-NAcount,NAcount/))
            R22=reshape(pack(estR,thisNA2),(/NAcount,NAcount/))
            invR11R12=R12
            call dposv('U',K-NAcount,NAcount,R11,K-NAcount,invR11R12,K-NAcount,info) !invR11R12=inv(R11)*R12
            
            !compute variance of imputed data
            call dgemm('T','N',NAcount,NAcount,K-NAcount,-1d0,invR11R12,K-NAcount,R12,K-NAcount,1d0,R22,NAcount) !R22=R22-invR11R12'*R12
         end if

         !compute means of imputed data
         E1=pack(Ey(i,:)-estE,.not.thisNA)
         E2=pack(estE,thisNA)
         call dgemm('T','N',NAcount,1,K-NAcount,1d0,invR11R12,K-NAcount,E1,K-NAcount,1d0,E2,NAcount) !E2=E2-invR11R12'*E1

         !assign means and variance of imputed data
         thisEy=Ey(i,:)
         thisEy=unpack(E2,thisNA,thisEy)
         thisEyy=0
         thisEyy=unpack(reshape(R22,(/NAcount**2/)),thisNA2,thisEyy)

      end if

      !add E^2[thisy] to V[thisy]
      call dgemm('N','T',K,K,1,1d0,thisEy,K,thisEy,K,1d0,thisEyy,K) !thisEyy=thisEyy+thisEy*thisEy'
      Ey(i,:)=thisEy
      Eyy=Eyy+thisEyy/N
      lastNA=thisNA

   end do iloop

   !update estimates (M-step)
   do j=1,K
      Enew(j)=sum(Ey(:,j))/N
   end do
   Rnew=Eyy
   call dgemm('N','T',K,K,1,-1d0,Enew,K,Enew,K,1d0,Rnew,K) !Rnew=Rnew-Enew*Enew'

   !check convergence
   Edist=maxval(abs(Enew-estE))
   Rdist=maxval(abs(Rnew-estR))
   write(*,'(A I4 A ES8.2)') 'iter ',iter,', dist=',max(Edist,Rdist)
   if ( Edist < tol .and. Rdist < tol ) then
      exit iterloop
   end if
   estE=Enew
   estR=Rnew

end do iterloop

!store results
if (iter>=maxiter) then
   write(*,'(A)') 'Reached maxiter'
else
   write(*,'(A)') 'Converged.'
end if

maxiter=iter
estE=Enew
estR=Rnew

end

