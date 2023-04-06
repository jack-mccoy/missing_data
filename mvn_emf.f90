! homemade EM to estimate mean and variance subject to missing data
! Fortran version

! in/outs 
!	X = data / imputed data (N by K)
! 	mu = guessed mu, estimated mu (K by 1)
!	Sig = guessed Sigma, estimated Sigma (K by K)
! ins
!	NAindex = indexes of missing values

! Algorithm
! E-step: 
!     For each obs i, compute 
!        e = E(X(i,:) | X(i,obs), mu, Sig)
!        S = E(X(i,:)*X(i,:)' | X(i,obs), mu, Sig)
!     Update X with e's
!     Summarize S with meanS
! M-step: compute mu_new and Sig_new
!     mu_new = mean(X)
!     Sig_new = meanS + mu_new*mu_new'

subroutine mvn_emf(X, mu, Sig, NAindex, tol, maxiter, K, N, update_mu)
	
implicit none

!declare input arguments	
integer K,N
double precision X(N,K), mu(K), Sig(K,K)
logical NAindex(N,K)
double precision tol
integer maxiter
logical update_mu ! .true. updates mu in each iteration (default).  .false. uses inputted mu throughout.

!declare other stuff
integer i,j,iter,NAcount,info
double precision mu_new(K), Sig_new(K,K)
double precision meanS(K,K), e(K), S(K,K)
double precision Sig_dist, mu_dist
logical thisNA(K), lastNA(K), thisNA2(K,K), thisDataNA(K,K), thisData2(K,K)
double precision, allocatable :: eo_less_muo(:), em(:)
double precision, allocatable :: Sig_oo(:,:), Sig_om(:,:), Vmm(:,:), invVooVom(:,:)

!calculations start here
write(*,'(A)') 'Starting EM:'

!allocate first, because we always deallocate first in the loop
allocate(eo_less_muo(1),em(1),Sig_oo(1,1),Sig_om(1,1),Vmm(1,1),invVooVom(1,1)) 

iterloop: do iter=1,maxiter

   meanS=0
   lastNA=.false.
   
   
   ! ==== compute means and variances for every observation (E-step) ====
   iloop: do i=1,N

      NAcount=count(NAindex(i,:))
      if (NAcount == 0) then 
         ! === complete data ===
         e=X(i,:)

         ! S := Xi*Xi'
         call dgemm('N','T',K,K,1,1d0,e,K,e,K,1d0,S,K)          

      elseif (NAcount == K) then          
         ! === no data ===
         e=mu

         ! S := Sig + mu*mu'
         S=Sig
         call dgemm('N','T',K,K,1,1d0,e,K,e,K,1d0,S,K) 

      else 
         ! === some data ===

         thisNA=NAindex(i,:)
                  
         if (any(thisNA .neqv. lastNA)) then 
            ! find "priors"
            ! can skip if missingness pattern is the same as last stock
			
            ! set up indexes   
            do j=1,K
               thisNA2(j,:)=thisNA(j) .and. thisNA
               thisDataNA(j,:)= (.not.thisNA(j)) .and. thisNA
               thisData2(j,:)= (.not.thisNA(j)) .and. (.not.thisNA)
            end do

            ! re-allocate
            deallocate(eo_less_muo,em,Sig_oo,Sig_om,Vmm,invVooVom)
            allocate(eo_less_muo(K-NAcount),em(NAcount))
            allocate(Sig_oo(K-NAcount,K-NAcount))
            allocate(Sig_om(K-NAcount,NAcount))
            allocate(Vmm(NAcount,NAcount))
            allocate(invVooVom(K-NAcount,NAcount))
            
            ! grab Sig components
            Sig_oo=reshape(pack(Sig,thisData2),(/K-NAcount,K-NAcount/))
            Sig_om=reshape(pack(Sig,thisDataNA),(/K-NAcount,NAcount/))            
            
            ! compute slope			
            ! solve for invVooVom in Sig_oo*invVooVom = Sig_om 
            ! arg 6 in dposv is in/out: in = Sig_om, out = invVooVom 
            invVooVom=Sig_om !arg 6 in dposv
            call dposv('U',K-NAcount,NAcount,Sig_oo,K-NAcount,invVooVom,K-NAcount,info) 
            
            ! compute Vmm = V(Xim|Xiobs)
            ! Vmm := Sig_mm - invVooVom'*Sig_om 
            ! netlib notation: C := A**T*B + C		 
   			!    dgemm(TRA,TRB,M,     ,N       ,K       ,ALPHA,A       ,LDA      ,B  ,K        ,BETA,C ,LDC)
            Vmm = reshape(pack(Sig,thisNA2),(/NAcount,NAcount/))            
            call dgemm('T','N',NAcount,NAcount,K-NAcount,-1d0,invVooVom,K-NAcount,Sig_om,K-NAcount,1d0,Vmm,NAcount) 
         end if		   

         ! replace missing 1st moments with expectations
         ! em :=  mum + invVooVom'*eo_less_muo 
		   ! netlib notation: C := C + A**T*B
		   !    dgemm(TRA,TRB,M,     ,N,K        ,ALPHA,A        ,LDA      ,B ,K        ,BETA,C ,LDC)
         em=pack(mu,thisNA)
         eo_less_muo=pack(X(i,:)-mu,.not.thisNA) ! 
         call dgemm('T','N',NAcount,1,K-NAcount,1d0  ,invVooVom,K-NAcount,eo_less_muo,K-NAcount,1d0 ,em,NAcount) 

		   !	e (1 x K) combines observed data in X(i,:) with em
         e=X(i,:)
         e=unpack(em,thisNA,e)
		 
		   ! initialize second moment S with Vmm
         S=0
         S=unpack(reshape(Vmm,(/NAcount**2/)),thisNA2,S)

         ! finish second moment S := S+e*e'
         call dgemm('N','T',K,K,1,1d0,e,K,e,K,1d0,S,K) 

      end if
     
	  
	  ! update
      X(i,:) = e ! store full imputed dataset
      meanS=meanS+S/N ! reduce second moments to just what we need
      lastNA=thisNA

   end do iloop

   ! ==== update estimates (M-step) ====
   
   ! Update mean 
   if ( update_mu ) then
      do j=1,K
         mu_new(j)=sum(X(:,j))/N
	   end do
   else
      mu_new = mu
   end if
   
   ! Update covariance
   ! Sig_new := meanS - mu_new*mu_new'
   ! Netlib Notation:
   ! C :=beta*C + alf*oA(A)*oB(B) 
   !          oA  oB  M N K  alf A    DA B  DB beta C   DC
   Sig_new=meanS   
   call dgemm('N','T',K,K,1,-1d0,mu_new,K,mu_new,K,1d0,Sig_new,K) 

   !check convergence
   mu_dist=maxval(abs(mu_new-mu))
   Sig_dist=maxval(abs(Sig_new-Sig))
   write(*,'(A I4 A ES8.2)') 'iter ',iter,', dist=',max(mu_dist,Sig_dist)
   if ( mu_dist < tol .and. Sig_dist < tol ) then
      exit iterloop
   end if
   mu=mu_new
   Sig=Sig_new

end do iterloop

!store results
if (iter>=maxiter) then
   write(*,'(A)') 'Reached maxiter'
else
   write(*,'(A)') 'Converged.'
end if

maxiter=iter
mu=mu_new
Sig=Sig_new

end

