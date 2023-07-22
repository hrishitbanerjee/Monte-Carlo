!                                                  
!3D Ising model - Metropolis Algorithm - Importance sampling     
!                                                  
! Monte Carlo code by Hrishit Banerjee


      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)

      parameter(idim1 = 25)

    dimension s(0:idim1, 0:idim1, 0:idim1),np(0:idim1, 0:idim1, 0:idim1)


      open(1, file='data.in', status='unknown')
      read(1,*) n, kount, kount_equil
       write(*,*)'min,max,istp'
       read(*,*)beg,fin,stp
       write(*,*)'J1'
       read(*,*) c1



      open(20,file='avg_eng.out')
      open(21,file='specific_heat.out')
      open(22,file='avg_mag.out')
      open(23,file='suscep.out')


      do t = beg,fin,stp

         print*, 'temp:', t

!     Initializations

         eng = 0.0
         omag = 0.0

!     initializing spins to random up or down state


         do i = 1, n
            do j = 1, n
               do k = 1, n

                  s(i,j,k) = 4 ! ordered state
                  
!                  call random_number(r)  ! disordered state
!                  if(r.le.0.5)then
!                     s(i,j,k) = 1
!                  else
!                     s(i,j,k) = -1
!                  endif
                  
               enddo
            enddo
         enddo

!     Initial energy and magnetization


         do i = 1, n
            do j = 1, n
               do k = 1, n
         
                  i1 = i + 1
                  j1 = j + 1
                  k1 = k + 1

                  if(i1.gt.n) i1 = i1 - n
                  if(j1.gt.n) j1 = j1 - n
                  if(k1.gt.n) k1 = k1 - n
               
                  eng = eng - (s(i,j,k)*c1*(s(i1,j,k)+s(i,j1,k)))
                  omag = omag + s(i,j,k)

               enddo
            enddo
         enddo

!     Flipping of spin and calculation of change in energy 

         ei = eng
         sum_e = 0
         s_sum_e = 0
         sum_m = 0
         s_sum_m = 0
      
         do loop = 1, kount
         
!     print*, 'loop:', loop, kount

            do ii = 1, n
               do jj = 1, n
                  do kk = 1, n

!     Choosing a spin randomly

                     call random_number(rani)
                     call random_number(ranj)
                     call random_number(rank)

                     i = n * rani + 1
                     j = n * ranj + 1
                     k = n * rank + 1
                  
                     s(i,j,k) = - s(i,j,k)
                     
                     i1 = i + 1
                     i0 = i - 1
                     
                     j1 = j + 1
                     j0 = j - 1
                     
                     k1 = k + 1
                     k0 = k - 1
                  
                     if(i1.gt.n) i1 = i1 - n
                     if(i0.eq.0) i0 = i0 + n

                     if(j1.gt.n) j1 = j1 - n
                     if(j0.eq.0) j0 = j0 + n

                     if(k1.gt.n) k1 = k1 - n
                     if(k0.eq.0) k0 = k0 + n

!     Calculation for change in energy

de = -2 * s(i,j,k) * c1* (s(i0, j, k) + s(i1,j,k) + s(i, j0, k) + s(i,j1,k) )

                     dm = 2 * s(i, j, k)


!     Importance Sampling

                     if(de.le.0) then
                        
                        ei = ei + de
                        omag = omag + dm
                        
                     else
                        
                        prob = exp(-de/t)
                        call random_number(r)

                        if(prob.gt.r)then
                           ei = ei + de
                           omag = omag + dm
                           
                        else
                           s(i, j, k) = - s(i, j, k)
                        endif
                        
                     endif
                     
                  enddo
               enddo
            enddo
            
            amag = abs(omag)
            
!     Calculation various physical quantiries after equilibration

            if(loop.gt.kount_equil)then
               sum_e = sum_e + ei
               s_sum_e = s_sum_e + (ei * ei)
               sum_m = sum_m + amag
               s_sum_m = s_sum_m + (omag *omag)
            endif
            
         enddo
         
         diff = kount - kount_equil

!     Averaging over ensemble
         
         a_sum_e = sum_e/diff
         a_s_sum_e = s_sum_e/diff
         a_sum_m = sum_m/diff
         a_s_sum_m = s_sum_m/diff
         
! 65      format(2(f18.10,1x))
         
!     Output Data 
         
         write(20,*) t, a_sum_e/(n*n*n)
         write(21,*) t,((a_s_sum_e)-(a_sum_e)**2)/(n*n*n*t*t)
         write(22,*) t, a_sum_m/(n*n*n)
         write(23,*) t, ((a_s_sum_m)-(a_sum_m)**2)/(n*n*n*t)
         
         
      enddo
      
      stop
      end
