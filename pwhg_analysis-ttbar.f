c Program to analyse tt~ production

c Initialise histograms

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer j,l
      character * 20 prefix
      integer lenocc
      external lenocc

		
		call inihists

      do j=1,3
         if(j.eq.1) then
            prefix='t'
         elseif(j.eq.2) then
            prefix='tb'
         elseif(j.eq.3) then
            prefix='ttbar'
C          elseif(j.eq.4) then
C             prefix='btop'
C          elseif(j.eq.5) then
C          	prefix='bbtop'
         endif
         l=lenocc(prefix)
         call bookupeqbins(prefix(1:l)//'_y'  ,0.04d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_eta',0.04d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_pt' ,2d0,0d0,500d0)
         call bookupeqbins(prefix(1:l)//'_m'  ,2d0,0d0,500d0)
      enddo
   	end

c Analysis subroutine
		subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_rad.h'
      include 'pwhg_weights.h'   ! KH added 17/8/16
      character * 6 whcprg      
      common/cwhcprg/whcprg
		data whcprg/'NLO   '/
      real * 8  dsig0            ! KH added 17/8/16
      real * 8  dsig(7)          ! KH added 17/8/16
      integer   nweights         ! KH added 17/8/16
      logical   iniwgts          ! KH added 17/8/16
      data      iniwgts/.true./  ! KH added 17/8/16
      save      iniwgts          ! KH added 17/8/16
      integer   ihep,jhep     	! HEPEVT index
      integer   id1,id2
c particle id numbers for identifying them in the event record
      integer 	 i_top,i_atop
c particle momenta
      real * 8  p_top(4),p_tb(4) 
c particle properties
      real * 8  y,eta,pt,mass

C - KH - 17/8/16 - added block from here ...
      if (iniwgts) then
         write(*,*) '*********************'
         if(whcprg.eq.'NLO') then
            write(*,*) ' NLO ANALYSIS      '
            weights_num=0
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) ' LHE ANALYSIS      '
         elseif(WHCPRG.eq.'HERWIG') then
            write(*,*) ' HERWIG ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA') then
            write(*,*) ' PYTHIA ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA8') then
            write(*,*) ' PYTHIA8 ANALYSIS   '
         endif
         write(*,*) '*********************'
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         iniwgts=.false.
      endif

      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
      endif
      if(sum(abs(dsig)).eq.0) return
C - KH - 17/8/16 - down to here ; copied from DYNNLOPS's pwhg_analysis-minlo.f

      i_top=0
      i_atop=0

c Find the tops from the event record

      do jhep=1,nhep
         if(idhep(jhep).eq.6) i_top = jhep
         if(idhep(jhep).eq.-6) i_atop = jhep
      enddo

      id1=idhep(1)
      id2=idhep(2)
      if(id1.eq.21) id1=0
      if(id2.eq.21) id2=0

      p_top=phep(1:4,i_top)
      p_tb=phep(1:4,i_atop)

      call getyetaptmass(p_top,y,eta,pt,mass)
      call filld('t_y',y,dsig)
      call filld('t_eta',eta,dsig)
      call filld('t_pt',pt,dsig)
      call filld('t_m',mass,dsig)

      call getyetaptmass(p_tb,y,eta,pt,mass)
      call filld('tb_y',y,dsig)
      call filld('tb_eta',eta,dsig)
      call filld('tb_pt',pt,dsig)
      call filld('tb_m',mass,dsig)

      call getyetaptmass(p_top+p_tb,y,eta,pt,mass)
      call filld('ttbar_y',y,dsig)
      call filld('ttbar_eta',eta,dsig)
      call filld('ttbar_pt',pt,dsig)
      call filld('ttbar_m',mass,dsig)

   	end

		subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end

      function phepDot(p_A,p_B)
      implicit none
      real * 8  phepDot
      real * 8  p_A(4),p_B(4)
      phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)
     1       -p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end

c KH - 17/8/16 - added the routines below as they are being called by
c the main-PYTHIA8.f code. They were added purely for this reason. I
c have no clue as to the physics behind doing it. More importantly,
c more generally, I have no idea what the main-PYTHIA8 and pythia8F77.cc
c are doing. 
      subroutine boost2reson4(pres,nm,pin,pout)
      implicit none
      integer nm
      real * 8 pres(4),pin(4,nm),pout(4,nm)
      real * 8 vec(3),beta
      beta=sqrt(pres(1)**2+pres(2)**2+pres(3)**2)/pres(4)
      vec(1)=pres(1)/(beta*pres(4))
      vec(2)=pres(2)/(beta*pres(4))
      vec(3)=pres(3)/(beta*pres(4))
      call mboost4(nm,vec,-beta,pin,pout)
      end


      
      subroutine mboost4(m,vec,beta,vin,vout)
c     boosts the m vectors vin(4,m) into the vectors vout(4,m) (that can
c     be the same) in the direction of vec(3) (|vec|=1) with velocity
c     beta.  Lorents convention: (t,x,y,z).
      implicit none
      integer m
      real * 8 vec(3),beta,vin(4,m),vout(4,m)
      real * 8 betav,gamma
      real * 8 vdotb
      integer ipart,idim
      gamma=1/sqrt(1-beta**2)
      do ipart=1,m
         vdotb=vin(1,ipart)*vec(1)
     #         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
         do idim=1,3
            vout(idim,ipart)=vin(idim,ipart)
     #           +vec(idim)*((gamma-1)*vdotb
     #           +gamma*beta*vin(4,ipart))
         enddo
         vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
      enddo
      end



