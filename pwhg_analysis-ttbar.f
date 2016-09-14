c Program to analyse tt~ production. Mainly copied/distilled 
c from pwhg_analysis-top-reconstruction.f in ttb_NLO_dec

c Initialise histograms

		subroutine init_hist
		implicit none
		include  'LesHouches.h'
		include 'pwhg_math.h'
		integer i,j,k,m,n,ls1,ls2,ls3,ls4,lp1
		character * 20 s1,s2,s3,s4 			! suffices for naming the histograms
		character * 20 p1 						! prefix for the histogram
		integer lenocc								! length of the character array
		external lenocc							! containing the prefix

		call inihists

		do i=1,2
			if(i.eq.1) s1 = '-incl'
			if(i.eq.2) s1 = '-str'
C 			if(i.eq.3) s1 = '-unstr'

			do j=1,3
				if(j.eq.1) s2 = '-all-angle'
				if(j.eq.2) s2 = '-coll'
				if(j.eq.3) s2 = '-wa'

				do k=1,2
					if(k.eq.1) s3 = '-all-mttb'
					if(k.eq.2) s3 = '-mttb-gt-1TeV'
C 					if(k.eq.2) s3 = '-mttb-gt-2TeV'

c (1.) Observables of the (i) top, (ii) t~, (iii) tt~ system
c (iv-vii) 1-4th hardest jets.

					do m=1,4
						if(m.eq.1) p1 = 't'
						if(m.eq.2) p1 = 'ttb'
						if(m.eq.3) p1 = 'j1'
						if(m.eq.4) p1 = 'j2'

						ls1=lenocc(s1)
						ls2=lenocc(s2)
						ls3=lenocc(s3)
						lp1=lenocc(p1)

						call bookupeqbins(p1(1:lp1)//'_y'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,-8d0,8d0)
						call bookupeqbins(p1(1:lp1)//'_eta'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,-8d0,8d0)
! different ranges on the pt spectra
						if(m.le.3) then
							call bookupeqbins(p1(1:lp1)//'_pt_2GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),2d0,0d0,500d0)
						else
							call bookupeqbins(p1(1:lp1)//'_pt_2GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),2d0,0d0,100d0)
						endif
						if(m.eq.3) then
							call bookupeqbins(p1(1:lp1)//'_pt_10GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),10d0,0d0,1000d0)
c							call bookupeqbins(p1(1:lp1)//'_pt_50GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),50d0,0d0,2000d0)
						endif
					enddo

c (2.) Rapidities in the y_ttbar=0 frame
	      		call bookupeqbins('yj1_minus_yttb'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,-8d0,8d0)
	      		call bookupeqbins('yj2_minus_yttb'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,-8d0,8d0)

c (3.) N additional jets
	      		call bookupeqbins('Njets_10GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),1d0,-0.5d0,10.5d0)
	      		call bookupeqbins('Njets_25GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),1d0,-0.5d0,10.5d0)
	      		call bookupeqbins('Njets_40GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),1d0,-0.5d0,10.5d0)

c (4.) distances between jets
c between the two light jets
	     			call bookupeqbins('dR_j1_j2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	     			call bookupeqbins('deta_j1_j2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	     			call bookupeqbins('dphi_j1_j2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),pi/50,0d0,pi)
cc between the two b jets
c	      		call bookupeqbins('dR_b1_b2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.08d0,0d0,8d0)
c	      		call bookupeqbins('deta_b1_b2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.08d0,0d0,8d0)
c	      		call bookupeqbins('dphi_b1_b2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),pi/50,0d0,pi)
c between the hardest b jet and the hardest light jet
	      		call bookupeqbins('dR_b1_j1'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	      		call bookupeqbins('deta_b1_j1'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	      		call bookupeqbins('dphi_b1_j1'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),pi/50,0d0,pi)
c between the hardest b jet and the lighter light jet
	      		call bookupeqbins('dR_b1_j2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	      		call bookupeqbins('deta_b1_j2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	      		call bookupeqbins('dphi_b1_j2'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),pi/50,0d0,pi)
c (5.) New observable that we have imagined
	      		call bookupeqbins('jet_pull-R'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	      		call bookupeqbins('jet_pull-eta'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),0.32d0,0d0,8d0)
	      		call bookupeqbins('jet_pull-phi'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),pi/50,0d0,pi)
c (6.) Total transverse momentum
	      		call bookupeqbins('j-Ht'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),50d0,0d0,2000d0)

c (7.) Gap fraction
c					call bookupeqbins('gap-fraction-5GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),5d0,0d0,500d0)
c					call bookupeqbins('gap-fraction-10GeV'//s1(1:ls1)//s2(1:ls2)//s3(1:ls3),10d0,0d0,1000d0)

	      	enddo
	      enddo
	   enddo

! for comparison with NLO plots

	   call bookupeqbins('j1_pt_2GeV-total',2d0,0d0,500d0)
	   call bookupeqbins('j1_pt_10GeV-total',10d0,0d0,1000d0)
	   call bookupeqbins('j1_pt_2GeV-gg',2d0,0d0,500d0)
	   call bookupeqbins('j1_pt_10GeV-gg',10d0,0d0,1000d0)
	   call bookupeqbins('j1_pt_2GeV-gg-str',2d0,0d0,500d0)
	   call bookupeqbins('j1_pt_10GeV-gg-str',10d0,0d0,1000d0)
	   call bookupeqbins('j1_pt_2GeV-gg-unstr',2d0,0d0,500d0)
	   call bookupeqbins('j1_pt_10GeV-gg-unstr',10d0,0d0,1000d0)

		end


c Analysis subroutine
		subroutine analysis(dsig0)
		implicit none
		include 'hepevt.h'
		include 'pwhg_math.h' 
		include 'LesHouches.h'
		include 'nlegborn.h'
		include 'pwhg_rad.h'
		include 'pwhg_weights.h'   	! KH added 17/8/16
		character * 7 whcprg      
		common/cwhcprg/whcprg
		data whcprg/'NLO    '/
		real * 8  dsig0            	! KH added 17/8/16
		real * 8  dsig(7)          	! KH added 17/8/16
		integer   nweights         	! KH added 17/8/16
		logical   iniwgts          	! KH added 17/8/16
		data      iniwgts/.true./  	! KH added 17/8/16
		save      iniwgts          	! KH added 17/8/16
		integer   ihep,jhep     		! HEPEVT index
		integer   id,id1,id2
		integer   ixx,jxx,kxx,lxx,j,mxx,mu 	! loop indices
! particle id numbers for identifying them in the event record
		integer 	 i_top,i_atop,i_bfromtop,i_abfromatop
		integer	 i_bjet,i_abjet,i_jet
		integer 	 bhadfromtop,bhadfromatop,counter,counter2,counter3,c4
		integer   i_b1,i_b2,i_j1,i_j2
		real * 8  njets10,njets25,njets40 		! For counting the n-additional jets
! particle momenta
		real * 8  p_top(4),p_tb(4),p_b(4),p_bb(4),p_hist(4),p_jet(4)
! particle observables
		real * 8  y,eta,pt,mass,mttbar,yttbar,ptttbar,ptjH,
     1          ptbhadfromtop,ptbhadfromatop,y_t,y_tb,deltay
     	logical   str,unstr
! distances between jets
		real * 8 dy,deta,dr,dphi
! Fastjet stuff
		integer   maxtracks,maxjets
		parameter (maxtracks=nmxhep,maxjets=20)
		integer mjets,jetvec(maxtracks),in_jet,
     1        ngenerations,ngen0,sonofid,binhadron
		logical   isForClustering(maxtracks),isbhadron
		external sonofid,binhadron,isbhadron,in_jet
		common/cngenerations/ngenerations
c Jet observables
		real * 8 j_kt(maxjets),j_eta(maxjets),j_rap(maxjets)
		real * 8 j_phi(maxjets),j_p(4,maxjets)
c names for the histograms
		character * 20 pf1,sf1,sf2,sf3 		! a prefix and suffices for the histograms
		logical cond1,cond2,cond3,cond4		! conditions for making the cuts for the histograms
		integer ls1,ls2,ls3,lp1					! length of the suffices and prefices
		integer lenocc
		external lenocc
		real * 8 bin_edge,binsize			! These are for making
		integer 	nbins,bxx,b2xx,cdan,c1,c2					! the gap fraction plot

		ngenerations = 4

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

c Initialise all numbers as 0
		i_top=0
		i_atop=0
		i_bfromtop=0
		i_abfromatop=0
		i_bjet=0
		i_abjet=0
		i_b1 = 0
		i_b2 = 0
		i_j1 = 0
		i_j2 = 0
		IsForClustering = .false.
		i_jet = 0
		p_jet(:) = 0
		p_hist(:) = 0

c Find the tops (and the bs that they decay into) 
c from the event record
		do jhep=1,nhep
			id=idhep(jhep)
			if(idhep(jhep).eq.6) i_top = jhep 		! find t
			if(idhep(jhep).eq.-6) i_atop = jhep 	! and t~
			id=abs(id)
			if(id.eq.5.or.id.eq.24) then
				if(min(sonofid(6,jhep),sonofid(-6,jhep)).lt.1d5) then
					if(idhep(jhep).eq.5) i_bfromtop = jhep			! b
					if(idhep(jhep).eq.-5) i_abfromatop = jhep		! b~
				endif
			endif
! Select jets in the NLO case so that we don't use Fastjet
         if(whcprg.eq.'NLO') then
            if(jhep.gt.2) then ! If it's a final state parton
               if(abs(idhep(jhep)).lt.6.or.idhep(jhep).eq.21) then
                  i_jet = jhep
               endif
            endif
         else
         	i_jet = 0
         endif
c for jets, using only final state particles excluding leptons
         if(isthep(jhep).eq.1.and.(id.lt.11.or.id.gt.16)) then
            IsForClustering(jhep) = .true.
         else
            IsForClustering(jhep) = .false.
         endif
      enddo

c Call Fastjet to build jets for PYTHIA case
      if(whcprg.ne.'NLO') then
      	mjets = maxjets
      	call buildjets(mjets,j_kt,j_eta,j_rap,j_phi,j_p,jetvec,
     1     	isForClustering)
      endif

      if(whcprg.ne.'NLO') then
c Find the b hadrons
     		bhadfromtop = 0
      	bhadfromatop = 0
      	ptbhadfromtop = 0
      	ptbhadfromatop = 0
c copied from ttb_NLO_dec
      	do j=1,nhep
         	if(IsForClustering(j).and.isbhadron(idhep(j))) then
            	if(binhadron(idhep(j)).eq.5) then ! is it a b or b~?
! Look for hardest (largest pt) hadron with a b quark content.
! Store in bhadfromtop, ptbhadfromtop
                	if(bhadfromtop.ne.0) then
c                 	 write(*,*) ' a top with more than one b son'
                  	call getyetaptmass(phep(1:4,j),y,eta,pt,mass)
                  	if(pt.gt.ptbhadfromtop) then
                     	bhadfromtop = j
                     	ptbhadfromtop = pt
            	      endif
               	else
                  	bhadfromtop = j
                  	call getyetaptmass(phep(1:4,j),y,eta,pt,mass)
                  	ptbhadfromtop = pt
               	endif
            	elseif(binhadron(idhep(j)).eq.-5) then
c same for bbar
      	         if(bhadfromatop.ne.0) then
c        	          write(*,*) ' a top with more than one b son'
            	      call getyetaptmass(phep(1:4,j),y,eta,pt,mass)
               	   if(pt.gt.ptbhadfromatop) then
                  	   bhadfromatop = j
                     	ptbhadfromatop = pt
            	      endif
               	else
                  	bhadfromatop = j
                  	call getyetaptmass(phep(1:4,j),y,eta,pt,mass)
                  	ptbhadfromatop = pt
               	endif
            	endif
         	endif
      	enddo

c Figure out which jets came from the b's i.e. which contain b hadrons
      	i_bjet = in_jet(bhadfromtop,jetvec)
      	i_abjet = in_jet(bhadfromatop,jetvec)
      	i_b1=min(i_bjet,i_abjet)
      	i_b2=max(i_bjet,i_abjet)
c We will require at least 2 of the 4 hardest jets to be b-jets (as in MC_TTBAR rivet analysis)
c If this is not the case, we shall skip the event
      	if(whcprg.eq.'LHE') then
      		counter=counter+1
      		if(i_b1.eq.1.and.i_b2.eq.2) then
      			i_j1=3
      		elseif(i_b1.eq.1.and.i_b2.eq.3) then
      			i_j1=2
      		elseif(i_b1.eq.2.and.i_b2.eq.3) then
      			i_j1=1
      		elseif(i_b1.eq.i_b2) then
      			counter2=counter2+1
      			write(25,*) 'Both bs clustered into the same jet'
	      		write(25,*) 'Skipping event:',counter2,'of',counter
	      		write(25,*)
	      		return
	      	endif
      	else
      		counter=counter+1
	      	if(i_bjet.gt.4.or.i_abjet.gt.4) then
	      		counter2=counter2+1
	      		write(25,*) 'We have less than two b-jets in our hardest 4 jets'
	      		write(25,*) 'Skipping event:',counter2,'of',counter
	      		write(25,*)
	      		return
	      	else
	      		if(i_bjet.eq.i_abjet) then
	      			counter2=counter2+1
	      			write(25,*) 'Both bs clustered into the same jet'
	      			write(25,*) 'Skipping event:',counter2,'of',counter
	      			write(25,*)
	      			return
	      		endif
	      		if(i_b1.eq.1.and.i_b2.eq.2) then
	      			i_j1=3
	      			i_j2=4
	      		elseif(i_b1.eq.1.and.i_b2.eq.3) then
	      			i_j1=2
	      			i_j2=4
	      		elseif(i_b1.eq.1.and.i_b2.eq.4) then
	      			i_j1=2
	      			i_j2=3
	      		elseif(i_b1.eq.2.and.i_b2.eq.3) then
	      			i_j1=1
	      			i_j2=4
	      		elseif(i_b1.eq.2.and.i_b2.eq.4) then
	      			i_j1=1
	      			i_j2=3
	      		elseif(i_b1.eq.3.and.i_b2.eq.4) then
	      			i_j1=1
	      			i_j2=2
	      		endif
	      	endif
	      endif
      endif

C       counter3=counter3+1
C  		write(23,*) 'event:', counter3
C       do ixx=1,20
C       	if(j_kt(ixx).gt.0) then
C       		if(i_bjet.eq.ixx) write(23,*) j_kt(ixx),'b jet'
C       		if(i_abjet.eq.ixx) write(23,*) j_kt(ixx),'b~ jet'
C       		if(i_bjet.ne.ixx.and.i_abjet.ne.ixx) write(23,*) j_kt(ixx)
C       	endif
C       enddo
C       write(23,*) 'b jet:',i_bjet
C       write(23,*) 'b~ jet:',i_abjet
C       write(23,*)

c Declare the momenta that will be used in the analysis
      id1=idhep(1)
      id2=idhep(2)
      if(id1.eq.21) id1=0
      if(id2.eq.21) id2=0
      p_top=phep(1:4,i_top)
      p_tb=phep(1:4,i_atop)
C       p_b=phep(1:4,i_bfromtop)
C       p_bb=phep(1:4,i_abfromatop)
      p_jet=phep(1:4,i_jet)

      njets10 = 0.0
      njets25 = 0.0
      njets40 = 0.0
      do j=1,mjets
      	if(j.ne.i_b1.and.j.ne.i_b2) then ! don't want to count b jets
         	if(j_kt(j).gt.10) then
            	njets10 = njets10 + 1.0
         	endif
         	if(j_kt(j).gt.25) then
            	njets25 = njets25 + 1.0
         	endif
         	if(j_kt(j).gt.40) then
            	njets40 = njets40 + 1.0
         	endif
         endif
      enddo

c calculate kinematic quantities needed for stretched/unstretched spectra
      call getyetaptmass(p_top,y,eta,pt,mass)
      y_t=y
      call getyetaptmass(p_tb,y,eta,pt,mass)
      y_tb=y
      deltay = y_t - y_tb
      call getyetaptmass(p_top+p_tb,y,eta,pt,mass)
      mttbar=mass
      yttbar=y
c Decide whether the event is a stretched or unstretched configuration
      str = .false.
      unstr = .false.
      if((rho.eq.1.and.deltay.lt.0).or.(rho.eq.2.and.deltay.gt.0)) then
      	str = .true.
      else
      	unstr = .true.
      endif

c Fill histograms - make the cuts

c First fill in the pT spectra for the NLO comparison

		if(whcprg.eq.'NLO') then
			p_hist(:)=p_jet(:)
		else
			p_hist(:)=j_p(:,i_j1)
		endif

		call getyetaptmass(p_hist,y,eta,pt,mass)

		if(pt.gt.0) then
			call filld('j1_pt_2GeV-total',pt,dsig)
			call filld('j1_pt_10GeV-total',pt,dsig)
			if(rho.eq.1.or.rho.eq.2) then
				call filld('j1_pt_2GeV-gg',pt,dsig)
				call filld('j1_pt_10GeV-gg',pt,dsig)
				if(str) then
					call filld('j1_pt_2GeV-gg-str',pt,dsig)
					call filld('j1_pt_10GeV-gg-str',pt,dsig)
				elseif(unstr) then
					call filld('j1_pt_2GeV-gg-unstr',pt,dsig)
					call filld('j1_pt_10GeV-gg-unstr',pt,dsig)
				endif
			endif
		endif


      do ixx=1,2

      	cond1 = .false.

      	if(ixx.eq.1) then
      		sf1 = '-incl'
      		cond1 = .true.
      	elseif(ixx.eq.2) then
      		sf1 = '-str'
      		if(str) then
      			cond1 = .true.
      		endif
      	endif

      	do jxx=1,3

      		cond2 = .false.

      		if(jxx.eq.1) then
      			sf2 = '-all-angle'
      			cond2 = .true.
      		elseif(jxx.eq.2) then
      			sf2 = '-coll'
      			if(coll.eq.1) then
      				cond2 = .true.
      			endif
      		elseif(jxx.eq.3) then
      			sf2 = '-wa'
      			if(coll.eq.0) then
      				cond2 = .true.
      			endif
      		endif

      		do kxx=1,2

      			cond3 = .false.

      			if(kxx.eq.1) then
      				sf3 = '-all-mttb'
      				cond3 = .true.
      			elseif(kxx.eq.2) then
      				sf3 = '-mttb-gt-1TeV'
      				if(mttbar.gt.1000) then
      					cond3 = .true.
      				endif
      			endif

c (1.) Observables of the (i) top, (ii) tt~ system
c (iii-iv) 1st two hardest non b-jets.

      			do mxx=1,4
      				cond4 = .false.
      				if(mxx.eq.1) then
      					pf1 = 't'
      					do mu=1,4
      						p_hist(mu) = p_top(mu)
      					enddo
      					cond4 = .true.
      				elseif(mxx.eq.2) then
      					pf1 = 'ttb'
      					do mu=1,4
      						p_hist(mu) = p_top(mu) + p_tb(mu)
      					enddo
      					cond4 = .true.
      				elseif(mxx.eq.3) then
      					pf1 = 'j1'
      					do mu=1,4
      						if(whcprg.eq.'NLO') then
      							p_hist(mu) = p_jet(mu)
      						else
	      						p_hist(mu) = j_p(mu,i_j1)
	      					endif
	      				enddo
	      				call getyetaptmass(p_hist,y,eta,pt,mass)
	      				if(pt.gt.0) then
	      					cond4 = .true.
	      				endif
	      			elseif(mxx.eq.4) then
	      				pf1 = 'j2'
	      				do mu=1,4
	      					p_hist(mu) = j_p(mu,i_j2)
	      				enddo
	      				call getyetaptmass(p_hist,y,eta,pt,mass)
	      				if(pt.gt.0) then
	      					cond4 = .true.
	      				endif
	      			endif

						ls1=lenocc(sf1)
						ls2=lenocc(sf2)
						ls3=lenocc(sf3)
						lp1=lenocc(pf1)

						if(cond1.and.cond2.and.cond3.and.cond4) then
							call getyetaptmass(p_hist,y,eta,pt,mass)
							call filld(pf1(1:lp1)//'_y'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),y,dsig)
				         call filld(pf1(1:lp1)//'_eta'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),eta,dsig)
		 	   	      call filld(pf1(1:lp1)//'_pt_2GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),pt,dsig)
		 	   	      if(mxx.eq.3) then
		 	   	      	call filld(pf1(1:lp1)//'_pt_10GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),pt,dsig)
C c (1.a) Gap fraction
C 								binsize=5d0
C 								nbins=100
C 								do bxx=1,nbins
C 									bin_edge=binsize*bxx
C 									if(pt.lt.bin_edge) then
C 										call filld('gap-fraction-5GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),binsize*bxx - 0.00001,binsize)
C 									endif
C 								enddo		
c		 	   	      	call filld(pf1(1:lp1)//'_pt_50GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),pt,dsig)
		 	   	      endif
		 	   	   endif
		 	   	enddo

c (2.) Rapidities in the y_ttbar=0 frame
		 	   	if(cond1.and.cond2.and.cond3) then
		 	   		call filld('yj1_minus_yttb'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),j_rap(i_j1)-yttbar,dsig)
		 	   		call filld('yj2_minus_yttb'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),j_rap(i_j2)-yttbar,dsig)

c (3.) N additional jets
		 	   		call filld('Njets_10GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),njets10,dsig)
      				call filld('Njets_25GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),njets25,dsig)
      				call filld('Njets_40GeV'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),njets40,dsig)

c (4.) distances between jets

c between the two light jets
      				call getdydetadphidr(j_p(:,i_j1),j_p(:,i_j2),dy,deta,dphi,dr)
      				call filld('dR_j1_j2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dr,dsig)
      				call filld('deta_j1_j2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),deta,dsig)
      				call filld('dphi_j1_j2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dphi,dsig)
cc between the two b jets
c      				call getdydetadphidr(j_p(:,i_b1),j_p(:,i_b2),dy,deta,dphi,dr)
c      				call filld('dR_b1_b2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dr,dsig)
c      				call filld('deta_b1_b2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),deta,dsig)
c      				call filld('dphi_b1_b2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dphi,dsig)
c between the b1 and j1
      				call getdydetadphidr(j_p(:,i_b1),j_p(:,i_j1),dy,deta,dphi,dr)
      				call filld('dR_b1_j1'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dr,dsig)
      				call filld('deta_b1_j1'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),deta,dsig)
      				call filld('dphi_b1_j1'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dphi,dsig)   
c between the b1 and j2
      				call getdydetadphidr(j_p(:,i_b1),j_p(:,i_j2),dy,deta,dphi,dr)
      				call filld('dR_b1_j2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dr,dsig)
      				call filld('deta_b1_j2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),deta,dsig)
      				call filld('dphi_b1_j2'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dphi,dsig)      			

c (5.) New observable that we have imagined
      				if(rho.eq.1.or.rho.eq.2) then ! only look at gg terms
      					if(abs(deltay).lt.0.2) then ! selects terms where Bfact ~ 0.5 i.e. the two matrix elements are roughly equal magnitudes
		      				if(( j_rap(i_j1) - y_t ).lt.0) then ! look in the lefthand hemisphere only
      							if(j_kt(i_j1).gt.10) then ! make a pT cut on the jet
      								call getdydetadphidr(j_p(:,i_j1),p_top,dy,deta,dphi,dr)
	      							call filld('jet_pull-R'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dr,dsig)
			      					call filld('jet_pull-eta'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),deta,dsig)
	      							call filld('jet_pull-phi'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),dphi,dsig)
	      						endif
	      					endif
	      				endif
	      			endif

c (6.) Total transverse momentum
						ptjH = j_kt(1)+j_kt(2)+j_kt(3)+j_kt(4) ! the sum of the 4 hardest jets kT
	      			call filld('j-Ht'//sf1(1:ls1)//sf2(1:ls2)//sf3(1:ls3),ptjH,dsig)

	      		endif
	      	enddo
	      enddo
	   enddo
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

		subroutine buildjets(mjets,kt,eta,rap,phi,pjet,jetvechep,
     1                                               isForClustering)
c     arrays to reconstruct jets
      implicit  none
      include  'hepevt.h'
      integer   maxtracks,maxjets
C       parameter (maxtracks=nmxhep,maxjets=nmxhep) ! <- think this is a mistake here 
      parameter (maxtracks=nmxhep,maxjets=20)
      integer   mjets,jetvechep(maxtracks)
      real * 8  kt(maxjets),eta(maxjets),rap(maxjets),
     1     phi(maxjets),pjet(4,maxjets)
      logical   isForClustering(maxtracks)
      real * 8  ptrack(4,maxtracks),pj(4,maxjets)
      integer   jetvec(maxtracks),itrackhep(maxtracks)
      integer   ntracks,njets,tmp1(maxtracks)
      integer   j,k,mu,ixx
      real * 8  r,palg,ptmin,pp,tmp
      integer sonofid
      external sonofid
C - Initialize arrays and counters for output jets
      ntracks=0
      kt=0
      eta=0
      rap=0
      phi=0
      pjet=0
      do ixx=1,maxtracks
      	ptrack(:,ixx)=0
      	jetvec(ixx)=0
      	pj(:,ixx)=0
      	jetvechep(ixx)=0
      	tmp1(ixx)=0
      enddo
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
      r     = 0.4d0
      ptmin = 0d0
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,jetvec)
      mjets=min(mjets,njets)
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
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
     1         +vin(2,ipart)*vec(2)+vin(3,ipart)*vec(3)
			do idim=1,3
				vout(idim,ipart)=vin(idim,ipart)
     1           +vec(idim)*((gamma-1)*vdotb
     2           +gamma*beta*vin(4,ipart))
			enddo
			vout(4,ipart)=gamma*(vin(4,ipart)+vdotb*beta)
		enddo
		end

		subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end



