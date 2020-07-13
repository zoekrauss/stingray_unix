#include <fintrf.h>
      SUBROUTINE stingray_src(u,t,iprec,inSetS,
     &      ghead,knx,kny,knodes,
     &      arcList,arcHead,kmx,kmy,kmz,kmaxfs,
     &      zhang,tf_line_integ,
     &      tf_anisotropy,a_r,a_t,a_p)

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
      
      REAL        pi
      PARAMETER ( pi       =3.14159265358979)

!    input variables

      INTEGER      knx,kny,kmx,kmy,kmz,kmaxfs,knodes   !
      REAL         u(knodes),t(knodes)                 !Graph1
      INTEGER      iprec(knodes)                       !Graph1
      INTEGER      arcList(kmaxfs,3)                   !ArcSet
      INTEGER      arcHead(4)
      REAL         ghead(8) 
      REAL         zhang(knx,kny)                       !Zdatum
      LOGICAL      tf_line_integ, tf_anisotropy          !Contrl
      LOGICAL      inSetS(knodes)                       !Intrvl
      REAL         a_r(knodes), a_t(knodes), a_p(knodes)      !Graph7

!     variables previously in common
    
      INTEGER      nx,ny,nz,nodes,nxy                  !Graph2
      REAL         g_origin_x,g_origin_y               !Graph2
      REAL         gx,gy,gz                            !Graph2
      INTEGER      minQ,maxQ                           !Intrvl
      REAL         dtmin,maxI                          !Intrvl
      INTEGER      mx,my,mz,nfs                        !ArcSet
      REAL         dzg(-kmz:kmz)                       !ArcSet
      REAL         sepxy2(-kmx:kmx,-kmy:kmy)           !ArcSet
      REAL         theta_fs(-kmx:kmx,-kmy:kmy,-kmz:kmz)  !ArcSet
      REAL         phi_fs(-kmx:kmx,-kmy:kmy)             !ArcSet

!    local variables

      INTEGER      iw,jw,kw,nw
      INTEGER      iv,jv,kv
      INTEGER      indx,numQ,indQ
      INTEGER      ifs
      INTEGER      lowmen(knodes)
      REAL         dt(kmx*kmy*kmz*6),dz,xxx
      INTEGER      nv(kmx*kmy*kmz*6)
      INTEGER      i,j,k
      REAL         dx,dy,dxy

!  Mex communication

      integer      icnt
      character*50 xstring

!   Declaration of mex variables (communication back to matlab)

      integer      n,m
      mwPointer    mexPrintf
      mwPointer    mexEvalString

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
