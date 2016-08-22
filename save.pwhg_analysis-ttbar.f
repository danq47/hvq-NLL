c Program to analyse tt~ production

c Initialise histograms

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer j,l1,l2,l3
      character * 20 prefix1,prefix2,prefix3 	! prefixes for naming the histograms
      integer lenocc										! length of the character array
      external lenocc									! containing the prefix

		
		call inihists

		call bookupeqbins('pt-t',2d0,0d0,300d0)

C c (1.) Observables of the (i) top, (ii) t~, (iii) tt~ system
C c (iv-vii) 1-4th hardest jets.

C       do j=1,7
C          if(j.eq.1) then
C             prefix1='t'
C          elseif(j.eq.2) then
C             prefix1='tb'
C          elseif(j.eq.3) then
C             prefix1='ttbar'
C          elseif(j.eq.4) then
C             prefix1='j1'
C          elseif(j.eq.5) then
C          	prefix1='j2'
C         	elseif(j.eq.6) then
C          	prefix1='j3'
C         	elseif(j.eq.7) then
C          	prefix1='j4'
C          endif
C          l1=lenocc(prefix1)
C          call bookupeqbins(prefix1(1:l1)//'_y'  ,0.04d0,-4d0,4d0)
C          call bookupeqbins(prefix1(1:l1)//'_eta',0.04d0,-4d0,4d0)
C          if(j.le.4) then
C          	call bookupeqbins(prefix1(1:l1)//'_pt_2GeV',2d0,0d0,500d0)
C          	call bookupeqbins(prefix1(1:l1)//'_pt_5GeV',5d0,0d0,500d0)
C          	call bookupeqbins(prefix1(1:l1)//'_pt_10GeV',10d0,0d0,1000d0)
C          	call bookupeqbins(prefix1(1:l1)//'_pt_50GeV',50d0,0d0,1000d0)
C          else
C          	call bookupeqbins(prefix1(1:l1)//'_pt_2GeV',2d0,0d0,250d0)
C          	call bookupeqbins(prefix1(1:l1)//'_pt_5GeV',5d0,0d0,250d0)
C          	call bookupeqbins(prefix1(1:l1)//'_pt_10GeV',10d0,0d0,500d0)
C          	call bookupeqbins(prefix1(1:l1)//'_pt_50GeV',50d0,0d0,500d0)
C          endif
C          call bookupeqbins(prefix1(1:l1)//'_mass',5d0,150d0,750d0)
C       enddo

C c (2.) Rapidities in the y_ttbar=0 frame
C       call bookupeqbins('yj1-minus-yttb',0.08d0,-8d0,8d0)
C       call bookupeqbins('yj2-minus-yttb',0.08d0,-8d0,8d0)
C       call bookupeqbins('yj3-minus-yttb',0.08d0,-8d0,8d0)
C       call bookupeqbins('yj4-minus-yttb',0.08d0,-8d0,8d0)

C c (3.) recoil of the jet off the ttbar system. In LHEF, this should be zero by definition
C       call bookupeqbins('pt-j1-minus-pt-ttbar',1d0,-100d0,100d0)

C c (4.) Gap fraction


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
      character * 7 whcprg      
      common/cwhcprg/whcprg
		data whcprg/'NLO    '/
      real * 8  dsig0            ! KH added 17/8/16
      real * 8  dsig(7)          ! KH added 17/8/16
      integer   nweights         ! KH added 17/8/16
      logical   iniwgts          ! KH added 17/8/16
      data      iniwgts/.true./  ! KH added 17/8/16
      save      iniwgts          ! KH added 17/8/16
      integer   ihep,jhep     	! HEPEVT index
      integer   id,id1,id2
      integer   ixx,jxx,kxx,lxx 			! loop indices
! particle id numbers for identifying them in the event record
      integer 	 i_top,i_atop,i_bfromtop,i_abfromatop,i_j1,i_j2,
     1 			 i_bjet,i_abjet,tmp(10),jet_position(10),i_j3,i_j4
     	integer 	 bhadfromtop,bhadfromatop
c particle momenta
      real * 8  p_top(4),p_tb(4),p_b(4),p_bb(4)
c particle observables
      real * 8  y,eta,pt,mass,mttbar,yttbar,ptttbar
c Fastjet stuff
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer mjets,jetvec(maxtracks)
      integer in_jet,ngenerations
      external in_jet
C       logical sonofid
      integer sonofid
      external sonofid
      logical   isForClustering(maxtracks)
      common/cngenerations/ngenerations
c Jet observables
      real * 8 j_kt(maxjets),j_eta(maxjets),j_rap(maxjets),
     1     j_phi(maxjets),j_p(4,maxjets)
C c conditions for analysis
C       logical condition1,condition2,condition3
c names for the histograms
      character * 20 prefix1,prefix2,prefix3
      integer l1,l2,l3
      real*8  p_hist(4),x
      integer lenocc
      external lenocc

		ngenerations = 4 	! parameter used in sonofhep

C - KH - 17/8/16 - added block from here ...
      if (iniwgts) then
         write(*,*) '*********************'
         if(whcprg.eq.'NLO    ') then
            write(*,*) ' NLO ANALYSIS      '
            weights_num=0
         elseif(WHCPRG.eq.'LHE    ') then
            write(*,*) ' LHE ANALYSIS      '
         elseif(WHCPRG.eq.'HERWIG ') then
            write(*,*) ' HERWIG ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA ') then
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
      i_bfromtop=0
      i_abfromatop=0
      i_bjet=0
      i_abjet=0
      bhadfromtop=0
      bhadfromatop=0
C       jet_position(1:10)=0
C       i_j1=0
C       i_j2=0
C       i_j3=0
C       i_j4=0
C       IsForClustering = .false.

c Find the particles from the event record

      do jhep=1,nhep
      	id=idhep(jhep)
         if(idhep(jhep).eq.6) i_top = jhep 		! find t
         if(idhep(jhep).eq.-6) i_atop = jhep 	! and t~
         id=abs(id)
         if(id.eq.5) then ! find b and b~ from t and t~ decay
            if(min(sonofid(6,jhep),sonofid(-6,jhep)).lt.1d5) then
            	if(idhep(jhep).eq.5) i_bfromtop = jhep 		! b
               if(idhep(jhep).eq.-5) i_abfromatop = jhep 	! b~
            endif
         endif
c for jets, using only final state particles excluding leptons
         if(isthep(jhep).eq.1.and.(id.lt.11.or.id.gt.16)) then
            IsForClustering(jhep) = .true.
         else
     	      IsForClustering(jhep) = .false.
         endif
      enddo

      id1=idhep(1)
      id2=idhep(2)
      if(id1.eq.21) id1=0
      if(id2.eq.21) id2=0

      p_top=phep(1:4,i_top)
      p_tb=phep(1:4,i_atop)
      if(i_bfromtop.ne.0) p_b=phep(1:4,i_bfromtop)
      if(i_abfromatop.ne.0) p_bb=phep(1:4,i_abfromatop)

c Call Fastjet to build jets
      mjets = maxjets
      call buildjets(mjets,j_kt,j_eta,j_rap,j_phi,j_p,jetvec,
     1     isForClustering)


      i_bjet = in_jet(i_bfromtop,jetvec) 		! Which jets come from the bs
      i_abjet = in_jet(i_abfromatop,jetvec)

c Calculate the hardest jet that did not originate from a b quark
c this should the showered version of the POWHEG hardest emission
      do ixx=1,6
      	tmp(ixx)=ixx 			! temporary array
      	jet_position(ixx)=0
      	if(i_bjet.eq.ixx.or.i_abjet.eq.ixx) then
      		tmp(ixx)=0			! delete the position of the b jets
      	endif
      enddo

      jxx=0 						! jet_position(n) now contains
      do ixx=1,4		! the n-th hardest non-b jet
 100		continue
      	if(tmp(ixx+jxx).ne.0) then
      		jet_position(ixx)=ixx+jxx
      	else
      		jxx=jxx+1
      		goto 100
      	endif
      enddo

      i_j1=jet_position(1)
      i_j2=jet_position(2)
      i_j3=jet_position(3)
      i_j4=jet_position(4)

c Fill histograms - make the cuts


      call getyetaptmass(p_top,y,eta,pt,mass)
      call filld('pt-t',pt,dsig)
C       do kxx=1,7

C       	p_hist=0

C          if(kxx.eq.1) then
C             prefix1 = 't'
C             p_hist = p_top
C          elseif(kxx.eq.2) then
C             prefix1 = 'tb'
C             p_hist = p_tb
C          elseif(kxx.eq.3) then
C             prefix1 = 'ttbar'
C             p_hist = p_top + p_tb
C          elseif(kxx.eq.4) then
C          	prefix1 = 'j1'
C           	p_hist = j_p(1:4,i_j1)
C          	p_hist=0
C          elseif(kxx.eq.5) then
C          	prefix1 = 'j2'
C          	p_hist = j_p(1:4,i_j2)
C          elseif(kxx.eq.6) then
C          	prefix1 = 'j3'
C          	p_hist = j_p(1:4,i_j3)
C          elseif(kxx.eq.7) then
C          	prefix1 = 'j4'
C          	p_hist = j_p(1:4,i_j4)
C          endif

C          l1=lenocc(prefix1)

C          call getyetaptmass(p_hist,y,eta,pt,mass)
C          if(kxx.eq.3) then
C          	mttbar=mass
C          	yttbar=y
C          	ptttbar=pt
C          endif
C          call filld(prefix1(1:l1)//'_y',y,dsig)
C          call filld(prefix1(1:l1)//'_eta',eta,dsig)
C          call filld(prefix1(1:l1)//'_pt_2GeV',pt,dsig)
C          call filld(prefix1(1:l1)//'_pt_5GeV',pt,dsig)
C          call filld(prefix1(1:l1)//'_pt_10GeV',pt,dsig)
C          call filld(prefix1(1:l1)//'_pt_50GeV',pt,dsig)
C          call filld(prefix1(1:l1)//'_mass',mass,dsig)

C       enddo

C !       x =j_rap(i_j1)
C !      all filld('j1_pt_2GeV',j_kt(i_j1),dsig)
C !		call bookupeqbins('yj1-minus-yttb',j_rap(i_j1)-yttbar,dsig)
C !		call bookupeqbins('yj2-minus-yttb',j_rap(i_j2)-yttbar,dsig)
C !		call bookupeqbins('yj3-minus-yttb',j_rap(i_j3)-yttbar,dsig)
C !		call bookupeqbins('yj4-minus-yttb',j_rap(i_j4)-yttbar,dsig)

C !		call bookupeqbins('pt-j1-minus-pt-ttbar',j_kt(i_j1)-ptttbar,dsig)








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

C       function sonofid(m,k)
C       implicit none
C       logical sonofid
C       integer m,k
C       include  'hepevt.h'
C       integer j,kcurr
C       integer ngenerations
C       common/cngenerations/ngenerations
C       kcurr=k
C       do j=1,ngenerations
C          if(abs(idhep(kcurr)).eq.m) then
C             sonofid = .true.
C             return
C          endif
C          kcurr = jmohep(1,kcurr)
C          if(kcurr.eq.0) then
C             sonofid = .false.
C             return
C          endif
C       enddo
C       sonofid=.false.
C       end

 		subroutine buildjets(mjets,kt,eta,rap,phi,pjet,jetvechep,
     1                                               isForClustering)
c     arrays to reconstruct jets
      implicit  none
      include  'hepevt.h'
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer   mjets,jetvechep(maxtracks)
      real * 8  kt(maxjets),eta(maxjets),rap(maxjets),
     1     phi(maxjets),pjet(4,maxjets)
      logical   isForClustering(maxtracks)
      real * 8  ptrack(4,maxtracks),pj(4,maxjets)
      integer   jetvec(maxtracks),itrackhep(maxtracks)
      integer   ntracks,njets
      integer   j,k,mu
      real * 8  r,palg,ptmin,pp,tmp
      logical sonofid
      external sonofid
C - Initialize arrays and counters for output jets
      ptrack = 0
      jetvec = 0
      ntracks=0
      pj = 0
      kt=0
      eta=0
      rap=0
      phi=0
      pjet=0
      jetvechep=0
C - Extract final state particles to feed to jet finder
      do j=1,nhep
         if(.not.isForClustering(j)) cycle
         if(ntracks.eq.maxtracks) then
            write(*,*) 'analyze: need to increase maxtracks!'
            write(*,*) 'ntracks: ',ntracks
            call exit(-1)
         endif
         ntracks=ntracks+1
         ptrack(:,ntracks) = phep(1:4,j)
         itrackhep(ntracks)=j
      enddo
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------- C
C - Run FastJet - C
C --------------- C
C - R = 0.7   radius parameter
C - f = 0.75  overlapping fraction
      palg  = -1
      r     = 0.7d0
      ptmin = 10d0
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,jetvec)
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=ptrack(mu,k)+pj(mu,jetvec(k))
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         kt(j)=sqrt(pjet(1,j)**2+pjet(2,j)**2)
         pp = sqrt(kt(j)**2+pjet(3,j)**2)
         eta(j)=0.5d0*log((pp+pjet(3,j))/(pp-pjet(3,j)))
         rap(j)=0.5d0*log((pjet(4,j)+pjet(3,j))/(pjet(4,j)-pjet(3,j)))
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      jetvechep = 0
      do j=1,ntracks
         jetvechep(itrackhep(j))=jetvec(j)
      enddo
      end

C       function in_jet(i_part,jetvec)
C       implicit none
C       include 'hepevt.h'
C       integer   maxtracks,maxjets
C       parameter (maxtracks=nmxhep,maxjets=20)
C       integer in_jet,jetvec(maxtracks),i_part
C       integer j
C       logical sonofhep
C       external sonofhep
C       do j=1,nhep
C          if(jetvec(j).ne.0) then
C             if(sonofhep(i_part,j)) then
C                in_jet = jetvec(j)
C                return
C             endif
C          endif
C       enddo
C       in_jet = 0
C       end

      function in_jet(i_part,jetvec)
      implicit none
      include 'hepevt.h'
      integer   maxtracks,maxjets
      parameter (maxtracks=nmxhep,maxjets=20)
      integer in_jet,jetvec(maxtracks),i_part
      integer j
      integer sonofhep
      external sonofhep
      do j=1,nhep
         if(jetvec(j).ne.0) then
            if(sonofhep(i_part,j).lt.1d5) then
               in_jet = jetvec(j)
               return
            endif
         endif
      enddo
      in_jet = 0
      end

c old
C       function sonofhep(m,k)
C       implicit none
C       logical sonofhep
C       integer m,k
C       include  'hepevt.h'
C       integer j,kcurr
C       integer ngenerations
C       common/cngenerations/ngenerations
C       kcurr=k
C       do j=1,ngenerations
C          if(kcurr.eq.m) then
C             sonofhep = .true.
C             return
C          endif
C          kcurr = jmohep(1,kcurr)
C          if(kcurr.eq.0) then
C             sonofhep = .false.
C             return
C          endif
C       enddo
C       sonofhep = .false.
C       end


      function sonofid(m,k)
c if k'th particle in hep common block
c is son of a particle with idhep=m returns
c how far back (along the jmohep sequence) the ancestor is.
c It looks back for no more than ngenerations levels.
c In case of failure it returns 1 000 000.
      implicit none
      integer, intent(in):: m,k
      integer sonofid,sonofid0
      sonofid = sonofid0(m,k,0)
      end

      recursive function sonofid0(m,k,level) result(result)
      implicit none
      integer, intent(in):: m,k,level
      integer result
      include  'hepevt.h'
      integer k1,k2,r1,r2
      integer ngenerations
      common/cngenerations/ngenerations
      if(level.gt.ngenerations) then
         result = 1000000
         return
      endif
      if(idhep(k).eq.m) then
         result = level
         return
      endif
      k1 = jmohep(1,k)
      k2 = jmohep(2,k)
      r1 = sonofid0(m,k1,level+1)
      r2 = sonofid0(m,k2,level+1)
      result = min(r1,r2)
      end

      function sonofhep(m,k)
      implicit none
      integer, intent(in):: m,k
      integer sonofhep,sonofhep0
      sonofhep = sonofhep0(m,k,0)
      end

      recursive function sonofhep0(m,k,level) result(result)
      implicit none
      integer, intent(in):: m,k,level
      integer result
      include  'hepevt.h'
      integer k1,k2,r1,r2
      integer ngenerations
      common/cngenerations/ngenerations
      if(level.gt.ngenerations) then
         result = 1000000
         return
      endif
      if(k.eq.m) then
         result = level
         return
      endif
      k1 = jmohep(1,k)
      k2 = jmohep(2,k)
      r1 = sonofhep0(m,k1,level+1)
      r2 = sonofhep0(m,k2,level+1)
      result = min(r1,r2)
      end

		function isbhadron(idhep)
      implicit none
      logical isbhadron
      integer idhep
      integer idigit
c all b hadrons have a 5 either a third (mesons) or fourth (barions) digit
      if(abs(idhep).eq.5.or.idigit(3,idhep).eq.5
     1     .or.idigit(4,idhep).eq.5) then
         isbhadron = .true.
      else
         isbhadron = .false.
      endif
      end

      function binhadron(idhep)
      implicit none
      integer binhadron
      integer idhep
      integer idigit
c all b hadrons have a 5 either a third (mesons) or fourth (barions) digit
      if(abs(idhep).eq.5) then
         binhadron = idhep
      elseif(idigit(4,idhep).eq.5) then
         binhadron = sign(5,idhep)
      elseif(idigit(4,idhep).eq.0.and.idigit(3,idhep).eq.5
c This line is to avoid to count b bbar resonances as having definite flavour
     1        .and.idigit(2,idhep).ne.5 ) then
         binhadron = - sign(5,idhep)
      else
         binhadron = 0
      endif
      end

      function idigit(ipos,inum)
      implicit none
      integer idigit,ipos,inum
      if(ipos.eq.0) then
         write(*,*) ' error: idigit(ipos.inum), ipos<=0!'
         call exit(-1)
      endif
      idigit = int(mod(abs(inum),10**ipos)/10**(ipos-1))
      end

