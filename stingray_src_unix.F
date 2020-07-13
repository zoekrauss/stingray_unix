      program stingray_src_unix


!      SUBROUTINE stingray_src(u,t,iprec,inSetS,
!     &      ghead,knx,kny,knodes,
!     &      arcList,arcHead,kmx,kmy,kmz,kmaxfs,
!     &      zhang,tf_line_integ,
!     &      tf_anisotropy,a_r,a_t,a_p)

! DIJKSTRA AND INTERVAL VERSION (see Klimes and Kvasnicka, 1994, G. J. Int.)
!    
! Calculates travel times and ray paths of first arrivals.
! Includes bathymetry/topography (see Toomey et al., 1994)
! Includes simple form of anisotropy (cos(2*theta)); 
! (see Dunn and Toomey, 1998 and Dunn et al. 2001).
!
! Modified from Rob Dunn's version of dkInt.f 
!
! Uses optimal size forward stars for quick calculations.
!
!    u:                Slowness field (s/km)
!
!    ghead:            Header for slowness model:
!                      g_origin_x = ghead(1):    x-coord of u(1,1,1), (km) 
!                      g_origin_y = ghead(2):    y-coord of u(1,1,1), (km) 
!                      (nx,ny,nz) = ghead(3:5):  Total number of nodes in (x,y,z)-directions
!                      (gx,gy,gz) = ghead(6:8):  Node spacing in (x,y,z)-directions, (km) 
!
!    inSetS:           Logical that is true when mininum time of node is known.
!
!    arcList:          List of indices in forward star
!    arcHead:          Header for arcList
!                      mx = arcHead(1)
!                      my = arcHead(2)
!                      mz = arcHead(3)
!                      nfs= arcHead(4)
!
!    zhang:            Bathymetry/topography, km.  
!                      z is positive upwards and relative to sea level.
!
!    tf_line_integ:    Logical, when true uses bresnham's algorithm.
!
!    tf_anisotropy:    Logical, when true do anisotropy 
!
!    t:                Defines shortest path travel times.
!    iprec:            Defines shortest paths by nodal number.
!
!    nw:               Node number under consideration
!    nv:               indices of FS(mx,my,mz,nw)
!

      implicit none

! sets dimensions of nnx_,nny_,nmx_,nmy_,nmz_,nmaxfs_,nnodes_

      include 'stingray_basedims.F' 
      
      REAL        pi
      PARAMETER ( pi       =3.14159265358979)

!    input variables

      INTEGER      knx,kny,kmx,kmy,kmz,kmaxfs,knodes   !
      REAL         u(nnodes_),t(nnodes_)                 !Graph1
      INTEGER      iprec(nnodes_)                       !Graph1
      INTEGER      arcList(nmaxfs_,3)                   !ArcSet
      INTEGER      arcHead(4)
      REAL         ghead(8) 
      REAL         zhang(nnx_,nny_)                       !Zdatum
      LOGICAL      tf_line_integ, tf_anisotropy          !Contrl
      LOGICAL      inSetS(nnodes_)                       !Intrvl
      REAL         a_r(nnodes_), a_t(nnodes_), a_p(nnodes_)      !Graph7
      REAL         dum, idum
      INTEGER      catch 

!     variables previously in common
    
      INTEGER      nx,ny,nz,nodes,nxy                  !Graph2
      REAL         g_origin_x,g_origin_y               !Graph2
      REAL         gx,gy,gz                            !Graph2
      INTEGER      minQ,maxQ                           !Intrvl
      REAL         dtmin,maxI                          !Intrvl
      INTEGER      mx,my,mz,nfs                        !ArcSet
      REAL         dzg(-nmz_:nmz_)                       !ArcSet
      REAL         sepxy2(-nmx_:nmx_,-nmy_:nmy_)           !ArcSet
      REAL         theta_fs(-nmx_:nmx_,-nmy_:nmy_,-nmz_:nmz_)  !ArcSet
      REAL         phi_fs(-nmx_:nmx_,-nmy_:nmy_)             !ArcSet

!    local variables

      INTEGER      iw,jw,kw,nw
      INTEGER      iv,jv,kv
      INTEGER      indx,numQ,indQ
      INTEGER      ifs
      INTEGER      lowmen(nnodes_)
      REAL         dt(nmx_*nmy_*nmz_*6),dz,xxx
      INTEGER      nv(nmx_*nmy_*nmz_*6)
      INTEGER      i,j,k
      REAL         dx,dy,dxy
!
!
!

! open ascii files containing the variables defined in stingray.m

! k7
      open( unit = 21, file = '.stingray_k7')
      read(21,*) dum
      knx = nint(dum)
      read(21,*) dum
      kny = nint(dum)
      read(21,*) dum
      kmx = nint(dum)
      read(21,*) dum
      kmy = nint(dum)
      read(21,*) dum
      kmz = nint(dum)
      read(21,*) dum
      kmaxfs = nint(dum)
      read(21,*) dum
      knodes = nint(dum)
      close(21)

       
      
      if ((knx>nnx_).or.(kny>nny_).or.(kmx>nmx_).or.(kmy>nmy_)) then
        catch = 0
        open (unit = 21, file = '.stingray_catch')
        write(21,*), catch
        close(21)
        stop
      end if
      if ((kmz > nmz_).or.(kmaxfs > nmaxfs_).or.(knodes >nnodes_ )) then
        catch = 0
        open (unit = 21, file = '.stingray_catch')
        write(21,*), catch
        close(21)
        stop
      end if


! iprec
      open( unit = 21, file = '.stingray_iprec')
      do i = 1, knodes
        read(21,*) dum
        iprec(i) = nint(dum)
      end do
      close(21)

! arcList
      open( unit = 21, file = '.stingray_arcList')
      do j = 1,3
        do i = 1,kmaxfs
          read(21,*) dum
          arcList(i,j) = nint(dum)
        end do
      end do
      close(21)

! arcHead
      open( unit = 21, file = '.stingray_arcHead')
      do i = 1,4
        read(21,*) dum
        arcHead(i) = nint(dum)
      end do
      close(21)

! u
      open( unit = 21, file = '.stingray_u')
      do i = 1,knodes
        read(21,*) u(i)
      end do
      close(21)


! t
      open( unit = 21, file = '.stingray_t')
      do i = 1,knodes
        read(21,*) t(i)
      end do
      close(21)

! ghead
      open( unit = 21, file = '.stingray_ghead')
      do i = 1,8
        read(21,*) ghead(i)
      end do
      close(21)

! zhang
      open( unit = 21, file = '.stingray_zhang')
      do j = 1,kny
        do i = 1,knx
          read(21,*) zhang(i,j)
        end do
      end do
      close(21)

! a_r
      open( unit = 21, file = '.stingray_a_r')
      do i = 1,knodes
        read(21,*) a_r(i)
      end do
      close(21)
      
      print*, 'right after reading in'
      print*, a_r(50)

! a_t
      open( unit = 21, file = '.stingray_a_t')
      do i = 1,knodes
        read(21,*) a_t(i)
      end do
      close(21)

! a_p
      open( unit = 21, file = '.stingray_a_p')
      do i = 1,knodes
        read(21,*) a_p(i)
      end do
      close(21)

! tf_line_integ
      open( unit = 21, file = '.stingray_tf_line_integ')
      read(21,*) dum
      idum = nint(dum)
      if (idum.eq.0) then
        tf_line_integ = .False.
      else
        tf_line_integ = .True.
      end if
      close(21)
     
      
      print*, 'tf_line_integ'
      print*,tf_line_integ
      

! tf_anisotropy
      open( unit = 21, file = '.stingray_tf_anisotropy')
      read(21,*) dum
      idum = nint(dum)
      if (idum.eq.0) then
        tf_anisotropy = .False.
      else
        tf_anisotropy = .True.
      end if
      close(21)

      
      print*, 'tf_anisotropy'
      print*,tf_anisotropy
      

! inSetS
      open( unit = 21, file = '.stingray_inSetS')
      do i = 1,knodes
        read(21,*) dum
        idum = nint(dum)
        if (idum.eq.0) then
          inSetS(i) = .False.
        else
          inSetS(i) = .True.
        end if
      end do
      close(21)

      
   

!
!
!



! Call to subroutine
      call stingray_src(u,t,iprec,inSetS,
     &                  ghead,knx,kny,knodes,
     &                  arcList,arcHead,kmx,kmy,kmz,kmaxfs,
     &                  zhang,tf_line_integ,
     &                  tf_anisotropy,a_r,a_t,a_p)

!
!
!

! Write output variables

! t
      open (unit = 21, file = '.stingray_t')
      do i = 1,knodes
        write(21,*) t(i)
      end do
      close(21)

! iprec
      open (unit = 21, file = '.stingray_iprec')
      do i = 1,knodes
        write(21,*) iprec(i)
      end do
      close(21)

! catch if everything went correctly
      catch = 1
      open (unit = 21, file = '.stingray_catch')
      write(21,*), catch
      close(21)

      end program






C----------------------------------------------------

      SUBROUTINE stingray_src(u,t,iprec,inSetS,ghead,knx,kny,
     &                        knodes,
     &                        arcList,arcHead,kmx,kmy,kmz,kmaxfs,
     &                        zhang,tf_line_integ,
     &                        tf_anisotropy,a_r,a_t,a_p)



      include 'stingray_basedims.F' ! sets dimensions of nnx_,nny_,nmx_,nmy_,nmz_,nmaxfs_,nnodes_


!    input variables

      INTEGER      knx,kny,kmx,kmy,kmz,kmaxfs,knodes   !
      REAL         u(nnodes_),t(nnodes_)                 !Graph1
      INTEGER      iprec(nnodes_)                       !Graph1
      INTEGER      arcList(nmaxfs_,3)                   !ArcSet
      INTEGER      arcHead(4)
      REAL         ghead(8) 
      REAL         zhang(nnx_,nny_)                       !Zdatum
      LOGICAL      tf_line_integ, tf_anisotropy          !Contrl
      LOGICAL      inSetS(nnodes_)                       !Intrvl
      REAL         a_r(nnodes_), a_t(nnodes_), a_p(nnodes_)      !Graph7
      REAL         dum, idum

!     variables previously in common
    
      INTEGER      nx,ny,nz,nodes,nxy                  !Graph2
      REAL         g_origin_x,g_origin_y               !Graph2
      REAL         gx,gy,gz                            !Graph2
      INTEGER      minQ,maxQ                           !Intrvl
      REAL         dtmin,maxI                          !Intrvl
      INTEGER      mx,my,mz,nfs                        !ArcSet
      REAL         dzg(-nmz_:nmz_)                       !ArcSet
      REAL         sepxy2(-nmx_:nmx_,-nmy_:nmy_)           !ArcSet
      REAL         theta_fs(-nmx_:nmx_,-nmy_:nmy_,-nmz_:nmz_)  !ArcSet
      REAL         phi_fs(-nmx_:nmx_,-nmy_:nmy_)             !ArcSet

!    local variables

      INTEGER      iw,jw,kw,nw
      INTEGER      iv,jv,kv
      INTEGER      indx,numQ,indQ
      INTEGER      ifs
      INTEGER      lowmen(nnodes_)
      REAL         dt(nmx_*nmy_*nmz_*6),dz,xxx
      INTEGER      nv(nmx_*nmy_*nmz_*6)
      INTEGER      i,j,k
      REAL         dx,dy,dxy

!  translate ghead values

      g_origin_x = ghead(1)
      g_origin_y = ghead(2)
      nx         = nint(ghead(3))
      ny         = nint(ghead(4))
      nz         = nint(ghead(5))
      gx         = ghead(6)
      gy         = ghead(7)
      gz         = ghead(8)
      nxy        = nx*ny
      nodes      = nx*ny*nz
    
C       n = mexPrintf ("Translated values")
C       n = mexPrintf ('\n')
C       n = mexPrintf ('\n')
C       n = mexEvalString('drawnow')
    
! Translate arcHead values

      mx         = arcHead(1)
      my         = arcHead(2)
      mz         = arcHead(3)
      nfs        = arcHead(4)
      
      

!  Note that theta and phi are not the same as hpt.
!  Initialize dzg and sepxy
!  Initialize theta & phi of forward star for 3D model 
!  Phi = 0 along x-axis.  Pos. theta is in +x,+y direction
!  Theta is measured from +z axis towards the x-y plane.
!  arcDist is faster and runs with regularly spaced grid: No topography

      DO k = -mz,mz,1

      !  Model is +ve upwards, graph k-level indices increase downward; flip sign (Oct 2009, DRT)
      !  Need to check that anisotropy is correct, given change in sign of dz/dzg.
         dz = -float(k)*gz
         dzg(k) = dz
         DO j = -my,my,1
            dy = float(j)*gy
            DO i = -mx,mx,1
               dx              = float(i)*gx
               dxy             = sqrt(dy**2 + dx**2)
               sepxy2(i,j)     = dy**2 + dx**2
               IF ((dx.EQ.0.).AND.(dy.EQ.0.)) THEN
                  phi_fs(i,j)   = 0.0
               ELSE
                  phi_fs(i,j)   = atan2(dy,dx)
               END IF
               IF ((dxy.EQ.0.).AND.(dz.EQ.0.)) THEN
                  theta_fs(i,j,k)   = 0.0
               ELSE
                  theta_fs(i,j,k)   = atan2(dxy,dz)
               END IF
               !arcDist(i,j,k)  = SQRT(dy**2 + dx**2 + dz**2)
            END DO
         END DO
      END DO
    
      ! Initialize dtmin
      ! Find minimum time arc

      dtmin = 10000.0
      DO i=1,nodes
         dtmin = min(u(i),dtmin)  
      END DO
      dtmin = AMIN1(gx,gy,gz)*dtmin
      maxI = 0.


    
      !!!!!!     Top of Main Loop     !!!!!!!

      minQ = 1
      maxQ = nodes

      ! icnt = 0

      DO WHILE (minQ .LE. maxQ)
         maxI = maxI + dtmin

         !  Update minQ and maxQ

         DO WHILE ((inSetS(minQ)) .AND. (minQ .LE. maxQ))
            minQ = minQ +1
         END DO
         DO WHILE ((inSetS(maxQ)) .AND. (maxQ .GE. minQ))
            maxQ = maxQ - 1
         END DO

         !  Selection (find tt in queue that are in lowest interval)

         numQ = 0
         DO indQ = minQ,maxQ
           IF (.NOT.inSetS(indQ) .AND. (t(indQ).LT.maxI)) THEN
              numQ = numQ + 1
              lowmen(numQ) = indQ
              inSetS(indQ) = .TRUE.
           END IF
         END DO

         DO indx = 1,numQ
            nw = lowmen(indx)
            kw  = ((nw-1)/nxy) + 1
            jw  = ((nw-(kw-1)*nxy-1)/nx) + 1
            iw  = nw - (jw-1)*nx - (kw-1)*nxy
            DO ifs = 1,nfs
               kv  = kw+arcList(ifs,3)
               jv  = jw+arcList(ifs,2)
               iv  = iw+arcList(ifs,1)
               IF ((kv .GT. 0) .AND. (kv .LE. nz) .AND.
     &             (jv .GT. 0) .AND. (jv .LE. ny) .AND.
     &             (iv .GT. 0) .AND. (iv .LE. nx)) THEN
                  nv(ifs) = (kv-1)*nxy+(jv-1)*nx+iv
                  IF (.NOT. inSetS(nv(ifs))) THEN
                     IF (tf_anisotropy) THEN

                       ! For anisotropy (can also do line integration)

                        CALL arcWeight(iw,jw,kw,iv,jv,kv,dt(ifs),
     &                                 u,a_r,a_t,a_p,ghead,knodes,
     &                                 zhang,knx,kny,
     &                                 sepxy2,dzg,theta_fs,phi_fs,
     &                                 kmx,kmy,kmz,tf_line_integ)

                     ELSE

                     ! No anisotropy/No line integration

                        dz = zhang(iv,jv) - zhang(iw,jw) + dzg(kv-kw)
                        xxx =sepxy2(arcList(ifs,1),arcList(ifs,2))+dz**2
                        dt(ifs) = (u(nw)+u(nv(ifs)))/2.0*SQRT(xxx)

                     END IF

                     IF (t(nw)+dt(ifs) .LT. t(nv(ifs))) THEN
                        iprec(nv(ifs)) = nw
                        t(nv(ifs))     = t(nw)+dt(ifs)
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO



      RETURN
      END


C----------------------------------------------------
       SUBROUTINE arcWeight(i,j,k,i2,j2,k2,dt,
     &                     u,a_r,a_t,a_p,ghead,knodes,
     &                     zhang,knx,kny,
     &                     sepxy2,dzg,theta_fs,phi_fs,
     &                     kmx,kmy,kmz,tf_line_integ)

!   arcWeight: rob dunn; modified drt
!
!   Usage:  To calculate weight between graph nodes
!
!   Description:
!
!   Bresenhams line algorithm, in 3 dimensions.
!       Includes 3D anisotropy.
!
!      (i,j,k)    - starting node
!      (i2,j2,k2) - ending node
!      dt         - travel time between the 2  nodes
!
!      First node assumed to be the origin
!
!   This function moves from node to node along the line between
!   two nodes 'n1' and 'n2'. Nodes on the line are selected using
!   Bresenhams algorithm).
!
!   See also:   Fundamentals of Interactive Computer Graphics
!           J.D. Foley and A. Van Dam
!           Addison Wesley Publishing Company,
!           Reading Mass.,
!           March 1983.
!           pp 433-436
!

      implicit none
      include 'stingray_basedims.F'
      
!  input variables

      INTEGER      knx,kny,kmx,kmy,kmz,knodes
      REAL         u(nnodes_),ghead(8)             
      REAL         a_r(nnodes_), a_t(nnodes_), a_p(nnodes_)
      REAL         zhang(nnx_,nny_)
      REAL         dzg(-nmz_:nmz_)
      REAL         sepxy2(-nmx_:nmx_,-nmy_:nmy_)
      REAL         theta_fs(-nmx_:nmx_,-nmy_:nmy_,-nmz_:nmz_)
      REAL         phi_fs(-nmx_:nmx_,-nmy_:nmy_)
      logical      tf_line_integ

!  variables previously in common
    
      INTEGER      nx,ny,nz,nodes,nxy                  !Graph2
      REAL         g_origin_x,g_origin_y               !Graph2
      REAL         gx,gy,gz                            !Graph2

!  local Variables

      INTEGER  i,j,k,i2,j2,k2
      INTEGER  inc,di,dj,dk,dx,dy,dz
      INTEGER  x,y,z,num,index
      INTEGER  x_step, y_step, z_step
      INTEGER  d1,d2
      INTEGER  incr1a,incr1b,incr2a,incr2b
      REAL     dt,real_inc
      REAL     sumu
      REAL     delz, dist,alpha

!  translate ghead values
      g_origin_x = ghead(1)
      g_origin_y = ghead(2)
      nx         = nint(ghead(3))
      ny         = nint(ghead(4))
      nz         = nint(ghead(5))
      gx         = ghead(6)
      gy         = ghead(7)
      gz         = ghead(8)
      nxy        = nx*ny
      nodes      = nx*ny*nz

!  Initialize slowness sum and node counter
!  Determine ijk step sizes

      sumu = 0.0
      inc  = 0

!  Note: always use di, dj, dk for calculating angles (this returns the
!  angle of the line segment of interest).
      di = i2-i
      dj = j2-j
      dk = k2-k
      dx = ABS(di)
      dy = ABS(dj)
      dz = ABS(dk)

!  Do the starting node

      num = (k-1)*nxy + (j-1)*nx + i
      CALL angle_spherical(a_t(num),a_p(num),
     &            theta_fs(di,dj,dk),phi_fs(di,dj),alpha)

!     &            theta_fs(di,dj),phi_fs(di,dj,dk),alpha)

      sumu  = u(num)/(1.+ (a_r(num)/2.0)*COS(2.0*alpha))

!  Do the intermediate node, provided line integration is true
      IF (((dx .GT. 2).OR.
     &     (dy .GT. 2).OR.
     &     (dz .GT. 2)).AND.tf_line_integ) THEN

         x = 0
         y = 0
         z = 0

         x_step = ISIGN(1,di)
         y_step = ISIGN(1,dj)
         z_step = ISIGN(1,dk)

!  Determine which direction is longest

         IF ((dx .GE. dy) .AND. (dx .GE. dz)) THEN

            incr1a = 2 * dy
            d1 = incr1a - dx
            incr1b = 2 * (dy - dx)

            incr2a = 2 * dz
            d2 = incr2a - dx
            incr2b = 2 * (dz - dx)

            DO index = 1,dx-1
               x = x + x_step
               IF (d1 .LT. 0) THEN
                  d1 = d1 + incr1a
               ELSE
                  y  = y  + y_step
                  d1 = d1 + incr1b
               END IF

               IF (d2 .LT. 0) THEN
                  d2 = d2 + incr2a
               ELSE
                  z  = z  + z_step
                  d2 = d2 + incr2b
               END IF

               inc = inc + 1
               num=(z+k-1)*nxy+(y+j-1)*nx+x+i
               CALL angle_spherical(a_t(num),a_p(num), 
     &                     theta_fs(di,dj,dk),phi_fs(di,dj),alpha)

!     &                      theta_fs(di,dj),phi_fs(di,dj,dk),alpha)
               sumu = sumu+(u(num)/(1.+(a_r(num)/2.0)*COS(2.0*alpha)))
            END DO

         ELSE IF ((dy .GT. dx) .AND. (dy .GT. dz)) THEN
            incr1a = 2 * dx
            d1 = incr1a - dy
            incr1b = 2 * (dx - dy)

            incr2a = 2 * dz
            d2 = incr2a - dy
            incr2b = 2 * (dz - dy)

            DO index = 1,dy-1
               y = y + y_step
               IF (d1 .LT. 0) THEN
                  d1 = d1 + incr1a
               ELSE
                  x  = x  + x_step
                  d1 = d1 + incr1b
               END IF
               IF (d2 .LT. 0) THEN
                  d2 = d2 + incr2a
               ELSE
                  z  = z  + z_step
                  d2 = d2 + incr2b
               END IF

               inc = inc + 1
               num=(z+k-1)*nxy+(y+j-1)*nx+x+i
               CALL angle_spherical(a_t(num),a_p(num), 
     &                     theta_fs(di,dj,dk),phi_fs(di,dj),alpha)

!     &                      theta_fs(di,dj),phi_fs(di,dj,dk),alpha)
               sumu = sumu+(u(num)/(1.+(a_r(num)/2.0)*COS(2.0*alpha)))
            END DO

! dz is the largest
         ELSE
            incr1a = 2 * dx
            d1 = incr1a - dz
            incr1b = 2 * (dx - dz)

            incr2a = 2 * dy
            d2 = incr2a - dz
            incr2b = 2 * (dy - dz)

            DO index = 1,dz-1
               z = z + z_step
               IF (d1 .LT. 0) THEN
                  d1 = d1 + incr1a
               ELSE
                  x  = x  + x_step
                  d1 = d1 + incr1b
               END IF
               IF (d2 .LT. 0) THEN
                  d2 = d2 + incr2a
               ELSE
                  y  = y  + y_step
                  d2 = d2 + incr2b
               END IF

               inc = inc + 1
               num=(z+k-1)*nxy+(y+j-1)*nx+x+i
               CALL angle_spherical(a_t(num),a_p(num), 
     &                     theta_fs(di,dj,dk),phi_fs(di,dj),alpha)
!     &                      theta_fs(di,dj),phi_fs(di,dj,dk),alpha)
               sumu = sumu+(u(num)/(1.+(a_r(num)/2.0)*COS(2.0*alpha)))
            END DO
         END IF

      ELSE
         inc = 0
      END IF

!  Do the last node
      num = (k2-1)*nxy + (j2-1)*nx + i2
      CALL angle_spherical(a_t(num),a_p(num), 
     &            theta_fs(di,dj,dk),phi_fs(di,dj),alpha)

!     &            theta_fs(di,dj),phi_fs(di,dj,dk),alpha)
      sumu = sumu + u(num)/(1.+ (a_r(num)/2.0) * COS(2.0*alpha))

!  Determine total number of nodes in the sum
      real_inc = FLOAT(2+inc)

!  Determine the arcWeight, dt

      delz   = zhang(i2,j2) - zhang(i,j) + dzg(k2-k)
      dist   = SQRT(sepxy2(i-i2,j-j2) + delz**2)

      dt = (sumu/real_inc)*dist

      RETURN
      END


C----------------------------------------------------
      SUBROUTINE angle_spherical(thet1,phi1,thet2,phi2,alpha)
!
!     Returns the angle between two vectors that are parameterized in
!     spherical coords. The length of the vectors is 1 in this routine.
!
!     theta is the zenith angle from the +z axis to the point (radians).
!
!     phi is the azimuth angle from the positive x-axis to the orthogonal
!     projection of the point in the x-y plane. Positive theta is in the +x,
!     +y direction.(degrees)
!
!     uses the definition of the dot product to find angle
!     11 trig, 9 arith
!
      REAL x1,x2,y1,y2,z1,z2
      REAL phi1,thet1,phi2,thet2,alpha


!   First convert to rectagular coords

      x1 = SIN(thet1)*COS(phi1)
      y1 = SIN(thet1)*SIN(phi1)
      z1 = COS(thet1)

      x2 = SIN(thet2)*COS(phi2)
      y2 = SIN(thet2)*SIN(phi2)
      z2 = COS(thet2)

      alpha = ACOS(x1*x2 + y1*y2 + z1*z2)

      RETURN
      END
