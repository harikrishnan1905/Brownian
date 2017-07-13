      implicit none
      integer i,istep,nstep,iseed
      double precision dt,gam,T,c0,c1,c2,sigr,sigv,crv,time,dx,dy,dz,
     .                 x,y,z,vx,vy,vz,zetarx,zetary,zetarz,zetavx,
     .                 zetavy,zetavz,g1,g2,Lx,Ly,Lz,
     .                 gasdev,ran2,x0,y0,z0,rsq(1000),fac
      external gasdev,ran2


      dt = 0.01d0
      gam = 100.0d0
      T = 1.0d0
      
      c0 = exp(-gam*dt)
      c1 = (1.0d0-c0)/(dt*gam)
      c2 = (1.0d0-c1)/(dt*gam)
      sigr = dt*T*(2.0d0-(3.0d0-4.0d0*exp(-gam*dt) 
     .     +exp(-2.0d0*gam*dt))/(dt*gam))/gam
      sigr = sqrt(sigr)
      sigv = T*(1.0d0-exp(-2.0d0*dt*gam))
      sigv = sqrt(sigv)
      crv = T*(1.0d0-exp(-gam*dt))
     .     *(1.0d0-exp(-gam*dt))/(gam*sigr*sigv)      

      iseed = -887
      !Lx = 10.0d0
      !Ly = 10.0d0
      !Lz = 10.0d0
      nstep = 1000
      vx = ran2(iseed)
      vy = ran2(iseed)
      vz = ran2(iseed)
      do i = 1,nstep
         rsq(i) = 0.0d0
      enddo

      do i = 1,1!0000

      print*,i 
      x0 = ran2(iseed)
      y0 = ran2(iseed)
      z0 = ran2(iseed)
      time = 0.0d0
      x = x0
      y = y0
      z = z0

      do istep = 1,nstep
        time = time + dt
        g1 = gasdev(iseed)
        g2 = gasdev(iseed)
        zetarx = sigr*g1
        zetavx = sigv*(crv*g1+sqrt(1.0d0-crv*crv)*g2)
        g1 = gasdev(iseed)
        g2 = gasdev(iseed)
        zetary = sigr*g1
        zetavy = sigv*(crv*g1+sqrt(1.0d0-crv*crv)*g2)
        g1 = gasdev(iseed)
        g2 = gasdev(iseed)
        zetarz = sigr*g1
        zetavz = sigv*(crv*g1+sqrt(1.0d0-crv*crv)*g2)
        x = x + c1*vx*dt + zetarx
        y = y + c1*vy*dt + zetary    
        z = z + c1*vz*dt + zetarz    
        write(100,*) x,y,z
        !x = x - Lx*anint((x/Lx) - 0.5d0)
        !y = y - Ly*anint((y/Ly) - 0.5d0)   
        !z = z - Lz*anint((z/Lz) - 0.5d0)   
        vx = c0*vx + zetavx
        vy = c0*vy + zetavy
        vz = c0*vz + zetavz
        dx = x0 - x
        dy = y0 - y
        dz = z0 - z
        !dx = dx - Lx*anint(dx/Lx)
        !dy = dy - Ly*anint(dy/Ly)
        !dz = dz - Lz*anint(dz/Lz)
        rsq(istep) = rsq(istep) + dx*dx + dy*dy + dz*dz
      enddo

      enddo

      time = 0.0d0
      do i = 1,istep
       time = time + dt
       fac = 6.0d0*(gam*time - 1.0d0 + exp(-time*gam))/(gam*gam)
       write(55,*)time,rsq(i)/10000.0d0,fac
      enddo
        
      stop
      end
