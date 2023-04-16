! homemade EM to estimate mean and variance subject to missing data
! Fortran version

! in/outs 
!	Ey = data / imputed data
! 	estE = guessed mu, estimated mu
!	estR = guessed Sigma, estimated Sigma
! ins
!	NAindex = indexes of missing values

subroutine mvn_emf(Ey, estE, estR, NAindex, tol, maxiter, K, N, update_estE)
use ieee_arithmetic
   
implicit none


!declare input arguments	
integer K,N
double precision Ey(N,K), estE(K), estR(K,K)
logical NAindex(N,K)
double precision tol
integer maxiter
logical update_estE ! .true. updates estE in each iteration (default).  .false. uses inputted estE throughout.

!declare other stuff
integer i,j,iter,NAcount,info
double precision Enew(K), Rnew(K,K)
double precision Eyy(K,K), thisEy(K), thisEyy(K,K)
double precision Rdist, Edist
logical thisNA(K), lastNA(K), thisNA2(K,K), thisDataNA(K,K), thisData2(K,K)
double precision, allocatable :: E1(:), E2(:)
double precision, allocatable :: R11(:,:), R12(:,:), R22(:,:), invR11R12(:,:)

!calculations start here
open(1, file = 'mvn_emf.log', status = 'unknown')  
write(1,'(A)') 'Starting EM:'

!allocate first, because we always deallocate first in the loop
allocate(E1(1),E2(1),R11(1,1),R12(1,1),R22(1,1),invR11R12(1,1)) 

iterloop: do iter=1,maxiter

   Eyy=0
   lastNA=.false.
   
   
   ! ==== compute means and variances for every observation (E-step) ====
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
		 
			! == prep calculations ==
			! 1 suffix indicates observed
			! 2 suffix indicates missing
			! computes
			!	slope: invR11R12' = R12'*inv(R11)
			!	new variance: R22 := R22 + invR11R12'R12			
			
            !set up indices            
            do j=1,K
               thisNA2(j,:)=thisNA(j) .and. thisNA
               thisDataNA(j,:)= (.not.thisNA(j)) .and. thisNA
               thisData2(j,:)= (.not.thisNA(j)) .and. (.not.thisNA)
            end do

            !re-allocate
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
            
			!compute slope			
			!solve for invR11R12 in R11*invR11R12 = R12 
			!	arg 6 in dposv is in/out: in = R12, out = invR11R12 
			invR11R12=R12 !arg 6 in dposv
            call dposv('U',K-NAcount,NAcount,R11,K-NAcount,invR11R12,K-NAcount,info) 
            
            !compute new variance 
			!R22 := invR11R12'*R12 + R22
			!netlib notation: C := A**T*B + C		 
			!    dgemm(TRA,TRB,M,     ,N       ,K       ,ALPHA,A       ,LDA      ,B  ,K        ,BETA,C ,LDC)
            call dgemm('T','N',NAcount,NAcount,K-NAcount,-1d0,invR11R12,K-NAcount,R12,K-NAcount,1d0,R22,NAcount) 
         end if

		 ! == "impute" missing values ==		 
		 ! impute with new means E2 :=  E2 + invR11R12'*E1 = E2 +  R12'*inv(R11)*E1
		 ! i.e. Ey(i,thisNA) = estE(thisNA) + estR(thisDataNA)'*estR(thisData2)^-1*(Ey(i,.not.thisNA)-estE(.not.thisNA))
         E1=pack(Ey(i,:)-estE,.not.thisNA) ! 
         E2=pack(estE,thisNA)
		 ! netlib notation: C := A**T*B + C		 
		 !    dgemm(TRA,TRB,M,     ,N,K        ,ALPHA,A        ,LDA      ,B ,K        ,BETA,C ,LDC)
         call dgemm('T','N',NAcount,1,K-NAcount,1d0  ,invR11R12,K-NAcount,E1,K-NAcount,1d0 ,E2,NAcount) 

         ! == reshape obs and imputed values into full vector ====
		 !	thisEy (1 x K) combines observed data in Ey(i,:) with means of missing data E2
         thisEy=Ey(i,:)
         thisEy=unpack(E2,thisNA,thisEy)
		 
		 ! == begin reshaping full second moment contribution ==
		 ! first reshape the variance correction
		 ! thisEyy (K x K) reshapes var of missing signals R22 to fit full model w/ 0s elsewhere, for now
		 !		thisEyy will be added to thisEy*thisEy' shortly		 
         thisEyy=0
         thisEyy=unpack(reshape(R22,(/NAcount**2/)),thisNA2,thisEyy)

      end if

	  ! == finish full second moment contribution ==	  
      ! add imputed second moment to variance correction
	  ! (note there is no correction if there are no missing values in the if statements above)
	  ! thisEyy := thisEyy+thisEy*thisEy'
      call dgemm('N','T',K,K,1,1d0,thisEy,K,thisEy,K,1d0,thisEyy,K) 
	  
	  ! update
      Ey(i,:)=thisEy ! store full imputed dataset
      Eyy=Eyy+thisEyy/N ! reduce second moment (each stock contributes thisEyy/N)
      lastNA=thisNA

   end do iloop

   ! ==== update estimates (M-step) ====
   
   ! Update mean 
   if ( update_estE ) then
	   do j=1,K
		  Enew(j)=sum(Ey(:,j))/N
	   end do
   else
	   Enew = estE
   end if
   
   ! Update covariance
   !Rnew := Eyy + Enew*Enew'
   Rnew=Eyy
   call dgemm('N','T',K,K,1,-1d0,Enew,K,Enew,K,1d0,Rnew,K) 

   !check convergence
   Edist=maxval(abs(Enew-estE))
   Rdist=maxval(abs(Rnew-estR))
   write(1,'(A I4 A ES8.2)') 'iter ',iter,', dist=',max(Edist,Rdist)

   if (max(Edist,Rdist) > 1e10) then
      write(1,'(A)') 'ERROR, DIVERGING'
      Enew = ieee_value(Enew, ieee_quiet_nan)
      Rnew = ieee_value(Rnew, ieee_quiet_nan)
      exit iterloop
   end if

   if ( Edist < tol .and. Rdist < tol ) then
      exit iterloop
   end if
   estE=Enew
   estR=Rnew

end do iterloop

!store results
if (max(Edist,Rdist) > 1e10) then
   write(1,'(A)') 'ERROR: Diverged'
else if (iter>=maxiter) then
   write(1,'(A)') 'ERROR: Reached maxiter'
else
   write(1,'(A)') 'Converged.'
end if

maxiter=iter
estE=Enew
estR=Rnew

close(1)

end

