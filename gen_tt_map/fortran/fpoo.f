c	ELLIPTICITY CORRECTIONS - Polynomial interpolation
c
      subroutine ellip()
      real sc0
      real sc1
      real sc2
c SQUARE ROOT OF 3. TIMES .5
      real s3
c EPICENTAL CO-LATITUDE - RADIANS
      real ecolat
c EPICENTRAL DISTANCE - RADIANS
      real edist
c AZIMUTH FROM EPICENTER TO RECEIVER - RADIANS
      real az
c EPICENTRAL DEPTH (km)
      real edepth
c ELLIPTICITY CORRECTION -OUTPUT
      real tcor
      character*(*) phase
c TAU's
      real t0, t1, t2
c	CONSTANTS FOR POLYNOMIAL
      real t0con(8,10), t1con(8,10), t2con(8,10)
      integer ii, j
      real adepth
      data (t0con(1,j),j=1,10) /-0.01711,-1.7791,0.,0.,0.,-0.9630,-13.
     &2326,13.7390,0.,0./
      data (t0con(2,j),j=1,10) /-0.08291,-2.1455,2.4538,-0.7907,0.,2.
     &0258,-12.9357,2.1287,5.2668,-0.9229/
      data (t0con(3,j),j=1,10) /-1.5022,-0.0943,1.9655,-1.1661,.1393,3.
     &4920,-9.9051,-0.3875,5.3581,-0.0686/
      data (t0con(4,j),j=1,10) /2.9971,-2.9549,0.4082,0.,0.,28.1650,9.
     &2160,-17.9030,-5.2995,3.2029/
      data (t0con(5,j),j=1,10) /3.6775,-2.2221,0.,0.,0.,-1.3127,-6.2476,
     &1.6684,0.,0./
      data (t0con(6,j),j=1,10) /-10.6238,15.4993,-7.4840,1.0673,0.,3.
     &2763,-6.4596,-0.4923,0.,0./
      data (t0con(7,j),j=1,10) /-0.01332,-3.2777,-1.2243,7.5246,0.,-3.
     &4856,-10.3187,43.4834,-70.5341,-50.2287/
      data (t0con(8,j),j=1,10) /-0.07859,-4.0924,4.6116, -1.4760,0.,2.
     &9104,-17.8661, 4.6262,7.1486,-1.9154/
      data (t1con(1,j),j=1,10) /.0040,-0.7841,6.0441,-17.5535,0.,-0.
     &2549,2.0519,-19.0605,-37.8235,54.5110/
      data (t1con(2,j),j=1,10) /-.0048, .0839,-2.2705,2.4137,-0.5957,-2.
     &4241,-4.2792,1.9728,3.5644,-0.5285/
      data (t1con(3,j),j=1,10) /.0033,-1.3485,0.1735,1.1583,-0.4162,-0.
     &1096,0.2576,-0.5978,0.1888,0.1600/
      data (t1con(4,j),j=1,10) /2.6249,-.0025,-0.2086,-0.0184,0.,-1.
     &5077,0.9904,0.3513,0.,0./
      data (t1con(5,j),j=1,10) /3.4213,-0.9359,0.,0.,0.,0.,0.,0.,0.,0./
      data (t1con(6,j),j=1,10) /-8.0633,8.0238,-1.7407,0.,0.,0.,0.,0.,0.
     &,0./
      data (t1con(7,j),j=1,10) /0.0109,-1.2300,8.9145,-27.5847,0.,-0.
     &6951,5.6201,-33.0908,-83.8233,102.4333/
      data (t1con(8,j),j=1,10) /-0.0311,0.1896,-4.0694,4.2599,-1.0387,-
     &3.9368,-8.4379,2.6814,6.9535,-0.6086/
      data (t2con(1,j),j=1,10) /0.0107,0.0275,-0.6912,0.0347,0.1157,-0.
     &1836,0.,0.0296,0.,0./
      data (t2con(2,j),j=1,10) /0.0107,0.0275,-0.6912,0.0347,0.1157,-0.
     &1836,0.,0.0296,0.,0./
      data (t2con(3,j),j=1,10) /0.0005,-0.01231,-1.0156,0.4396,0.,0.,0.,
     &0.,0.,0./
      data (t2con(4,j),j=1,10) /-3.5838,2.1474,-0.3548,0.,0.,-1.3369,-5.
     &4889,0.6809,1.5096,-0.0763/
      data (t2con(5,j),j=1,10) /-2.9912,1.0313,0.,0.,0.,0.,0.,0.,0.,0./
      data (t2con(6,j),j=1,10) /3.2814,-7.1224,3.5418,-0.5115,0.,0.,0.,
     &0.,0.,0./
      data (t2con(7,j),j=1,10) /0.00025,0.1685,-2.2435,3.3433,0.,-0.
     &0503,0.5353,1.5362,-14.3118,-3.2938/
      data (t2con(8,j),j=1,10) /0.0843,-0.2917,-0.6767,-0.2934,0.2779,-
     &0.4336,0.0306,0.07113,0.,0./
      data s3 /.8660254/
      save sc0, sc1, sc2
c
c	INITIAL CALL TO SET UP CONSTANTS
      entry ellref(ecolat)
c
      sc0 = .25*(1.+3.*cos(2.*ecolat))
      sc1 = s3 * sin(2.*ecolat)
      sc2 = s3 * sin(ecolat)*sin(ecolat)
      return
c
c	CALLED ONCE FOR EACH STATION - RETURNS ELLIPTICITY CORRECTION IN tcor
      entry ellcor(edist, azim, edepth, phase,tcor)
      adepth = edepth/6371.
c DETERMINE INDEX FOR POLYNOMIAL
      if(.not.(phase .eq. 'P'))goto 23000
         if(.not.(edist .lt. (15.*(3.14159265/180.))))goto 23002
            ii = 1
            goto 23003
c        else
23002       continue
            ii = 2
23003    continue
         goto 23001
c     else
23000    continue
         if(.not.(phase .eq. 'PcP'))goto 23004
            ii = 3
            goto 23005
c        else
23004       continue
            if(.not.(phase .eq. 'PKPab'))goto 23006
               ii = 4
               goto 23007
c           else
23006          continue
               if(.not.(phase .eq. 'PKPbc'))goto 23008
                  ii = 5
                  goto 23009
c              else
23008             continue
                  if(.not.(phase .eq. 'PKIKP'))goto 23010
                     ii = 6
                     goto 23011
c                 else
23010                continue
                     if(.not.(phase .eq. 'S'))goto 23012
                        if(.not.(edist .lt. (15.*(3.14159265/180.))))
     &                    goto 23014
                           ii = 7
                           goto 23015
c                       else
23014                      continue
                           ii = 8
23015                   continue
                        goto 23013
c                    else
23012                   continue
                        tcor = 0.
                        return
23013                continue
23011             continue
23009          continue
23007       continue
23005    continue
23001 continue
c	COMPUTE TAU's
c      write(6,*) ii
      t0 = t0con(ii,1) + edist*(t0con(ii,2) + edist*(t0con(ii,3) + 
     &edist*(t0con(ii,4) + edist*t0con(ii,5)))) +adepth*(t0con(ii,6) + 
     &adepth*t0con(ii,7)) + adepth*edist*(t0con(ii,8) + t0con(ii,9)*
     &adepth + t0con(ii,10)*edist)
      t1 = t1con(ii,1) + edist*(t1con(ii,2) + edist*(t1con(ii,3) + 
     &edist*(t1con(ii,4) + edist*t1con(ii,5)))) +adepth*(t1con(ii,6) + 
     &adepth*t1con(ii,7)) + adepth*edist*(t1con(ii,8) + t1con(ii,9)*
     &adepth + t1con(ii,10)*edist)
      t2 = t2con(ii,1) + edist*(t2con(ii,2) + edist*(t2con(ii,3) + 
     &edist*(t2con(ii,4) + edist*t2con(ii,5)))) +adepth*(t2con(ii,6) + 
     &adepth*t2con(ii,7)) + adepth*edist*(t2con(ii,8) + t2con(ii,9)*
     &adepth + t2con(ii,10)*edist)
      tcor = sc0 * t0 + sc1 * cos(azim) * t1 + sc2 * cos(2.*azim) * t2
      return
      end


C ELLIPTICITY CORRECTIONS FOR AK135 MODEL (full set of phases)
C
       SUBROUTINE kellip()
C==========================================================================
C                                                                         
C    Ellipticity correction for any given phase using
C    Dziewonski & Gilbert representation
C                                                   
C      The ellipticity corrections are found by linear interpolation       
C    in terms of values calculated for the ak135 model for a wide 
C    range of phases to match the output of the iasp software 
C
C     first call:  ellref(ecolat) 
C                        - to set up source dependent constants
C     2nd call  :  ellcor(phase,edist,depth,ecolat,azim,tcor,abrt) 
C                        - to calculate correction for a station                                                                                     C                                                                         
C    Parameters: 
C    character  
C          phase : a  string specifying the PHASE,   -e.g P, ScP etc.  
C                                                        
C    real 
C          edist  :  epicentral distance to station (in degrees)     
C          edepth :  depth of event         
C          ecolat :  epicentral co-latitude of source (in radians) 
C          azim   :  azimuth from source to station (in radians)
C                                
C          tcor   :  time correction for path to allow for ellipticity
C 
C    logical 
C          abrt   :  a logical variable -usally set to .FALSE.  
C                    which is set to .TRUE. if a phase for      
C                    which no data is available is chosen       
C                                                                         
C==========================================================================
C   B.L.N. Kennett RSES,ANU        May 1995, August 1996                 
C   (based on earlier routine by D.J. Brown)
C   with input from W. Spakman, Utrecht
C=========================================================================
      character *(*) phase
      character*8 phcod(57)
      integer phind(57),phspn(57),phnch(57)
      real edist,edepth,ecolat,azim,
     ^     sc0,sc1,sc2,s3,tcor,
     ^     tau0, a0,b0,h0,d0,e0,f0,g0,
     ^     tau1, a1,b1,h1,d1,e1,f1,g1,
     ^     tau2, a2,b2,h2,d2,e2,f2,g2
      real dpth(6),delta(50)
      real t0(50,6),t1(50,6),t2(50,6)
      integer Ne,Nd
      logical abrt
      data phcod/
     & "Pup   ","P     ","Pdiff ","PKPab ","PKPbc ","PKPdf ",
     & "PKiKP ","pP    ","pPKPab","pPKPbc","pPKPdf","pPKiKP",
     & "sP    ","sPKPab","sPKPbc","sPKPdf","sPKiKP","PcP   ",
     & "ScP   ","SKPab ","SKPbc ","SKPdf ","SKiKP ","PKKPab",
     & "PKKPbc","PKKPdf","SKKPab","SKKPbc","SKKPdf","PP    ",
     & "P'P'  ","Sup   ","S     ","Sdiff ","SKSac ","SKSdf ",
     & "pS    ","pSKSac","pSKSdf","sS    ","sSKSac","sSKSdf",
     & "ScS   ","PcS   ","PKSab ","PKSbc ","PKSdf ","PKKSab",
     & "PKKSbc","PKKSdf","SKKSac","SKKSdf","SS    ","S'S'  ",
     & "SP    ","PS    ","PnS   "/
      data phind/
     &        1,      14,      91,     136,     165,     178,
     &      235,     364,     433,     462,     475,     532,
     &      661,     742,     771,     784,     841,     970,
     &     1047,    1100,    1113,    1134,    1195,    1316,
     &     1337,    1382,    1507,    1516,    1573,    1702,
     &     1827,    1932,    1945,    2022,    2067,    2132,
     &     2197,    2234,    2295,    2356,    2425,    2490,
     &     2551,    2628,    2681,    2694,    2711,    2772,
     &     2781,    2838,    2967,    3140,    3273,    3398,
     &     3587,    3656,    3697/
      data phspn/
     &        3,      19,      11,       7,       3,      14,
     &       32,      17,       7,       3,      14,      32,
     &       20,       7,       3,      14,      32,      19,
     &       13,       3,       5,      15,      30,       5,
     &       11,      31,       2,      14,      32,      31,
     &       26,       3,      19,      11,      16,      16,
     &        9,      15,      15,      17,      16,      15,
     &       19,      13,       3,       4,      15,       2,
     &       14,      32,      43,      33,      31,      47,
     &       17,      10,       6/ 
      data phnch/
     &        3,       1,       5,       5,       5,       5,
     &        5,       2,       6,       6,       6,       6,
     &        2,       6,       6,       6,       6,       3,
     &        3,       5,       5,       5,       5,       6,
     &        6,       6,       6,       6,       6,       2,
     &        4,       3,       1,       5,       5,       5,
     &        2,       6,       6,       2,       6,       6,
     &        3,       3,       5,       5,       5,       6,
     &        6,       6,       6,       6,       2,       4,
     &        2,       2,       3/ 
      data dpth/ 0.0, 100.0, 200.0, 300.0, 500.0, 700.0 /
      save sc0,sc1,sc2
c...
c     In addition to the phase names listed above a number of phase aliases
c     are available in the routine phase_alias, e.g. Pn --> P etc
c     The input phase code is first checked against the phcod array
c     and next against the phase aliases.
c<sc>
c	                     initial call to set up source dependent constants
      entry kellref(ecolat)
c                                            
      s3 = sqrt(3.0)/2.0
      sc0 = 0.25*(1.0+3.0*cos(2.0*ecolat))
      sc1 = s3*sin(2.0*ecolat)
      sc2 = s3*sin(ecolat)*sin(ecolat)
      return
c<sc>
c<ec>                                           phase identification
      entry kellcor(phase,edist,edepth,ecolat,azim,tcor,abrt)
*      write(6,*) "phase,edist,edepth,ecolat,azim"
*      write(6,*)  phase,edist,edepth,ecolat,azim
      Nd = 6
      NUMPH = 57
      deldst = 5.0
      abrt = .FALSE.
c                                             check on the length of phase
      l=len(phase)
      if(l.lt.8) then
       stop 
     >    'character variable `phase` should have at least length 8'
      endif

c                                             select phase
      ip = -1
      nc=min(lnblk(phase),8)
      do 10 i=1,NUMPH
        if(nc.ne.phnch(i)) goto 10
        if (phase(1:nc) .eq. phcod(i)(1:nc)) then
          ip = i
          go to 11
        endif
 10   continue
 11   continue

      if(ip.eq.-1) then
c                                             check phase aliases
        call phase_alias(phase,edist,ip)
      endif
      Ne = phspn(ip)
*      write(6,*) "ip:",ip
c                                              phase not found
      if(ip.lt.0) then
        write(6,*) phase,'  is not available'
        abrt = .true.
        return
      endif
c                                               special case of upgoing waves
*
c
c                                                acquire phase information
       nr = phind(ip)
*       write(6,*) "nrec:",nr
       read(21,61,rec=nr) phcod(ip),np,d1,d2
*       write(6,*) "phcode,np,d1,d2: ", phcod(ip),np,d1,d2
       nr = nr+1
       if(np.ne.Ne) write(6,*) "HELP! - index wrong"
       do 15 i=1,np
         read(21,62,rec=nr) delta(i)
         nr = nr+1
         read(21,63,rec=nr) (t0(i,m),m=1,6)
         nr = nr+1
         read(21,63,rec=nr) (t1(i,m),m=1,6)
         nr = nr+1
         read(21,63,rec=nr) (t2(i,m),m=1,6)
         nr = nr+1
 15    continue         
 61    format(a8,i10,2f10.0)
 62    format(f10.0)
 63    format(6f10.4)
c                                  distance index
       idist = 1 + int( (edist-d1)/ deldst )
       if(edist.lt.d1) idist =1
       if(edist.gt.d2) idist= np-1
c                                  depth index
       do 25 j = 1,Nd-1
         if ((dpth(j).le.edepth).and.(dpth(j+1).ge.edepth))then
            jdepth = j
            goto 26
         endif
 25    continue
 26    continue
*       write(6,*) "idist, jdepth;",idist,jdepth
c
*                      need to allow for zero entries (where phase
*                      description strongly depth dependent)
c tau0
         a0 = t0(idist,jdepth)
         b0 = t0(idist,jdepth+1)
         h0 = t0(idist+1,jdepth+1)
         d0 = t0(idist+1,jdepth)
         e0 = a0 + 
     ^       (d0-a0)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f0 = b0 + 
     ^       (h0-b0)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g0 = e0 + 
     ^       (f0-e0)*(edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau0 = g0
c tau1
         a1 = t1(idist,jdepth)
         b1 = t1(idist,jdepth+1)
         h1 = t1(idist+1,jdepth+1)
         d1 = t1(idist+1,jdepth)
         e1 = a1 + 
     ^       (d1-a1)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f1 = b1 + 
     ^       (h1-b1)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g1 = e1 + 
     ^       (f1-e1)*(edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau1 = g1
c tau2
         a2 = t2(idist,jdepth)
         b2 = t2(idist,jdepth+1)
         h2 = t2(idist+1,jdepth+1)
         d2 = t2(idist+1,jdepth)
         e2 = a2 + 
     ^       (d2-a2)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f2 = b2 + 
     ^       (h2-b2)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g2 = e2 + 
     ^       (f2-e2)*(edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau2 = g2
c
*         write(6,*) "tau0,tau1,tau2:",tau0,tau1,tau2
         caz = cos(azim)
         cbz = cos(2.0*azim)
*         write(6,*) "azim,caz,cbz",azim,caz,cbz    
c
         tcor = sc0*tau0 + sc1*cos(azim)*tau1 + sc2*cos(2.0*azim)*tau2
c
      return
c<ec>
      end
      subroutine phase_alias(phase,delta,ip)

c-    check for alternative phase names
c     input phase, delta
c     output ip (index of phcod)

      character*(*) phase
      if(phase(1:3).eq.'Pg ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sg ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPg ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPg ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSg ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSg ') then
c       phase='sS      '
        ip=40
c
      elseif(phase(1:3).eq.'Pb ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sb ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPb ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPb ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSb ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSb ') then
c       phase='sS      '
c
      elseif(phase(1:3).eq.'Pn ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sn ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPn ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPn ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSn ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSn ') then
c       phase='sS      '
        ip=40
      else if(phase(1:4).eq.'SPn ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SPb ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SPg ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SnP ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'PSn ') then
c       phase='PS      '
        ip=56
      else if(phase(1:5).eq.'PnPn ') then
c       phase='PP      '
        ip=30
      else if(phase(1:5).eq.'SnSn ') then
c       phase='SS      '
        ip=53
c                                       upgoing P, S
      else if(phase(1:2).eq.'p ') then
c       phase='Pup     '
        ip=1  
      else if(phase(1:2).eq.'s ') then
c       phase='Sup     '
        ip=32 
c                                        
      else if(delta.le.100.0.and.phase.eq.'pPdiff  ') then
c       phase='pP      '
        ip=8
      else if(delta.le.100.0.and.phase.eq.'sPdiff  ') then
c       phase='sP      '
        ip=13
      else if(delta.le.100.0.and.phase.eq.'pSdiff  ') then
c       phase='pS      '
        ip=37
      else if(delta.le.100.0.and.phase.eq.'sSdiff  ') then
c       phase='sS      '
        ip=40
      else if(delta.ge.100.0.and.phase.eq.'pPdiff  ') then
c       phase='Pdiff   '
        ip=8
      else if(delta.ge.100.0.and.phase.eq.'sPdiff  ') then
c       phase='Pdiff   '
        ip=13
      else if(delta.ge.100.0.and.phase.eq.'pSdiff  ') then
c       phase='Sdiff   '
        ip=37
      else if(delta.ge.100.0.and.phase.eq.'sSdiff  ') then
c       phase='Sdiff    '
        ip=40
c            
      else if(delta.le.165.0.and.phase.eq.'PKPdiff ') then
c       phase='PKPbc '
        ip=5
      else if(delta.le.165.0.and.phase.eq.'pPKPdiff') then
c       phase='pPKPbc '
        ip=10
      else if(delta.le.165.0.and.phase.eq.'sPKPdiff') then
c       phase='sPKPbc '
        ip=15
c                             
      else if(phase(1:4).eq."P'P'") then
c       phase="P'P'    "
        ip =31
      else if(phase(1:4).eq."S'S'") then
c       phase="P'P'    "
        ip =54
c                            diffractions (approx)
      else if(delta.gt.100.0.and.phase.eq.'pPdiff  ') then
c       phase='Pdiff   '
        ip=8
      else if(delta.gt.100.0.and.phase.eq.'sPdiff  ') then
c       phase='Pdiff   '
        ip=13
      else if(delta.gt.100.0.and.phase.eq.'pSdiff  ') then
c       phase='Sdiff   '
        ip=37
      else if(delta.gt.100.0.and.phase.eq.'sSdiff  ') then
c       phase='Sdiff    '
c
      else
        ip=-1
      endif
      return
      end
c-----
      integer function lnblk(s)
      character *(*) s
      l = len(s)
      do i=l,1,-1
        if (s(i:i).gt.' ') then
          lnblk=i
          return
        endif
      end do   
      lnblk=0
      return
      end
      subroutine warn(msg)
      character*(*) msg
      write(*,100) msg
 100  format(1x,a)
      return
      end
      subroutine tnoua(ia,nc)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine tnoua writes the first nc charcters of the
c   character string ia to the standard 
c   output without the trailing newline (allowing user input on the 
c   same line).  Programmed on 17 September 1980 by R. Buland.
c
      save
      character*(*) ia
      write(*,100)ia(1:nc)
 100  format(a,$)
      return
      end
      subroutine dasign(lu,mode,ia,len)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine dasign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by dasign are
c   assumed to be direct access.  Programmed on 3 December 1979 by
c   R. Buland.
c
      save
      character*(*) ia
      logical exst
c
      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.3) ns=3
      go to (1,2),nf
 1    go to (11,12,13),ns
 11   open(lu,file=ia,status='old',form='formatted',
     1 access='direct',recl=len)
      return
 12   inquire(file=ia,exist=exst)
      if(exst) go to 11
 13   open(lu,file=ia,status='new',form='formatted',
     1 access='direct',recl=len)
      return
 2    go to (21,22,23),ns
 21   open(lu,file=ia,status='old',form='unformatted',access='direct',
     1 recl=len)
      return
 22   inquire(file=ia,exist=exst)
      if(exst) go to 21
 23   open(lu,file=ia,status='new',form='unformatted',access='direct',
     1 recl=len)
      return
      end
      subroutine vexit(ierr)
c      call exit(ierr)
      if (ierr .NE.0) print *,'exit called with error code',ierr
      stop
      end
      subroutine tabin(in,modnam)
      include 'ttlim.inc'
      character*(*) modnam
c     logical log
      character*8 phcd,phdif(6)
      character cdum*20
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
c
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
c      data tauc,xc/jtsm*0d0,jxsm*0d0/
c

      nin=in
      phdif(1)='P'
      phdif(2)='S'
      phdif(3)='pP'
      phdif(4)='sP'
      phdif(5)='pS'
      phdif(6)='sS'
c++
      call asnag1(nin,-1,1,'Enter model name:',modnam)

c++ 
      read(nin) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km,fcs,nafl,
     1 indx,kndx
      read(nin) pm,zm,ndex
      read(nin) pu,pux
      read(nin) phcd,px,xt,jndx
      read(nin) pt,taut
      read(nin) coef
      call retrns(nin)

c
      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)
      cdum=modnam(1:nb)//'.tbl'

      call dasign(nin,-1,cdum,nasgr)
c
      do 11 nph=1,2
 11   pu(ku(nph)+1,nph)=pm(1,nph)
c
c     write(10,*)'nasgr nl len2',nasgr,nl,len2
c     write(10,*)'nseg nbrn mt ku km',nseg,nbrn,mt,ku,km
c     write(10,*)'xn pn tn',xn,pn,tn
c     write(10,200)(i,(ndex(i,j),pm(i,j),zm(i,j),j=1,2),i=1,mt(2))
c200  format(/(1x,i3,i7,2f12.6,i7,2f12.6))
c     write(10,201)(i,(pu(i,j),j=1,2),i=1,ku(2)+1)
c201  format(/(1x,i3,2f12.6))
c     write(10,201)(i,(pux(i,j),j=1,2),i=1,km(2))
c     write(10,202)(i,(nafl(i,j),j=1,3),(indx(i,j),j=1,2),(kndx(i,j),
c    1 j=1,2),(fcs(i,j),j=1,3),i=1,nseg)
c202  format(/(1x,i3,7i5,3f5.0))
c     cn=180./3.1415927
c     write(10,203)(i,(jndx(i,j),j=1,2),(px(i,j),j=1,2),(cn*xt(i,j),
c    1 j=1,2),phcd(i),i=1,nbrn)
c203  format(/(1x,i3,2i5,2f12.6,2f12.2,2x,a))
c     write(10,204)(i,pt(i),taut(i),(coef(j,i),j=1,5),i=1,jout)
c204  format(/(1x,i5,0p2f12.6,1p5d10.2))
c
      tn=1./tn
      dn=3.1415927/(180.*pn*xn)
      odep=-1.
      ki=0
      msrc(1)=0
      msrc(2)=0
      k=1
      do 3 i=1,nbrn
      jidx(i)=jndx(i,2)
      do 4 j=1,2
 4    dbrn(i,j)=-1d0
 8    if(jndx(i,2).le.indx(k,2)) go to 7
      k=k+1
      go to 8
 7    if(nafl(k,2).gt.0) go to 9
      ind=nafl(k,1)
      l=0
      do 10 j=jndx(i,1),jndx(i,2)
      l=l+1
 10   tp(l,ind)=pt(j)
 9    if(nafl(k,1).gt.0.and.(phcd(i)(1:1).eq.'P'.or.
     1 phcd(i)(1:1).eq.'S')) go to 3
      do 5 j=1,6
      if(phcd(i).eq.phdif(j)) go to 6
 5    continue
      go to 3
 6    dbrn(i,1)=1d0
      phdif(j)=' '
 3    continue
c     write(10,205)(i,phcd(i),(dbrn(i,j),j=1,2),jidx(i),i=1,nbrn)
c205  format(/(1x,i5,2x,a,2f8.2,i5))
c     write(10,206)(i,(tp(i,j),j=1,2),i=1,jbrnu)
c206  format(/(1x,i5,2f12.6))
      return
      end
      subroutine depset(dep,usrc)
      save 
      include 'ttlim.inc'
      logical dop,dos,segmsk,prnt
      character*8 phcd
      real usrc(2)
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/prtflc/segmsk(jseg),prnt(2)
c      data segmsk,prnt/jseg*.true.,2*.false./
c

      if(amax1(dep,.011).ne.odep) go to 1
      dop=.false.
      dos=.false.
      do 2 i=1,nseg
      if(.not.segmsk(i).or.iidx(i).gt.0) go to 2
      if(iabs(nafl(i,1)).le.1) dop=.true.
      if(iabs(nafl(i,1)).ge.2) dos=.true.
 2    continue
      if(.not.dop.and..not.dos) return
      go to 3
c
 1    nph0=0
      int0(1)=0
      int0(2)=0
      mbr1=nbrn+1
      mbr2=0
      dop=.false.
      dos=.false.
      do 4 i=1,nseg
      if(.not.segmsk(i)) go to 4
      if(iabs(nafl(i,1)).le.1) dop=.true.
      if(iabs(nafl(i,1)).ge.2) dos=.true.
 4    continue
      do 5 i=1,nseg
      if(nafl(i,2).gt.0.or.odep.lt.0.) go to 5
      ind=nafl(i,1)
      k=0
      do 15 j=indx(i,1),indx(i,2)
      k=k+1
 15   pt(j)=tp(k,ind)
 5    iidx(i)=-1
      do 6 i=1,nbrn
 6    jndx(i,2)=-1
      if(ki.le.0) go to 7
      do 8 i=1,ki
      j=kk(i)
 8    pt(j)=pk(i)
      ki=0
c   Sample the model at the source depth.
 7    odep=amax1(dep,.011)
      rdep=dep
      if(rdep.lt..011) rdep=0.
      zs=amin1(alog(amax1(1.-rdep*xn,1e-30)),0.)
      hn=1./(pn*(1.-rdep*xn))
      if(prnt(1).or.prnt(2)) write(10,100)dep
 100  format(/1x,'Depth =',f7.2/)
c

 3    if(nph0.gt.1) go to 12
      if(dop) call depcor(1)
      if(dos) call depcor(2)
      go to 14
 12   if(dos) call depcor(2)
      if(dop) call depcor(1)
c
c   Interpolate all tau branches.
c
 14   j=1
      do 9 i=1,nseg
      if(.not.segmsk(i)) go to 9
      nph=iabs(nafl(i,1))
c     print *,'i iidx nph msrc nafl =',i,iidx(i),nph,msrc(nph),nafl(i,1)
      if(iidx(i).gt.0.or.(msrc(nph).le.0.and.nafl(i,1).gt.0)) go to 9
      iidx(i)=1
      if(nafl(i,2).le.0) int=nafl(i,1)
      if(nafl(i,2).gt.0.and.nafl(i,2).eq.iabs(nafl(i,1)))
     1  int=nafl(i,2)+2
      if(nafl(i,2).gt.0.and.nafl(i,2).ne.iabs(nafl(i,1)))
     1  int=iabs(nafl(i,1))+4
      if(nafl(i,2).gt.0.and.nafl(i,2).ne.nafl(i,3)) int=nafl(i,2)+6
 11   if(jndx(j,1).ge.indx(i,1)) go to 10
      j=j+1
      go to 11
 10   idel(j,3)=nafl(i,1)
c      print *,'spfit:  j int =',j,int 
      call spfit(j,int)
      mbr1=min0(mbr1,j)
      mbr2=max0(mbr2,j)
      if(j.ge.nbrn) go to 9
      j=j+1
c     print *,'j jidx indx jndx =',j,jidx(j),indx(i,2),jndx(j,2)
      if(jidx(j).le.indx(i,2).and.jndx(j,2).gt.0) go to 10
 9    continue
c     write(10,*)'mbr1 mbr2',mbr1,mbr2
c     write(10,*)'msrc isrc odep zs us',msrc,isrc,odep,sngl(zs),
c    1 sngl(us(1)),sngl(us(2))
c     write(10,200)ki,(i,iidx(i),kk(i),pk(i),i=1,nseg)
c200  format(/10x,i5/(1x,3i5,f12.6))
      usrc(1)=us(1)/pn
      usrc(2)=us(2)/pn

      return
      end
c
c----------------------------------------------------------------
c
      subroutine depcor(nph)
      save
      include 'ttlim.inc'
      character*8 phcd
      logical noend,noext,segmsk,prnt
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp,ua,taua
      double precision tup(jrec),umod,zmod,tauus1(2),tauus2(2),xus1(2),
     1 xus2(2),ttau,tx,sgn,umin,dtol,u0,u1,z0,z1,fac,du
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/pdec/ua(5,2),taua(5,2),deplim,ka
      common/prtflc/segmsk(jseg),prnt(2)
      equivalence (tauc,tup)
c      data tol,dtol,deplim,ka,lpower/.01,1d-6,1.1,4,7/ 
      data tol,dtol,lpower/.01,1d-6,7/
c
c     write(10,*)'depcor:  nph nph0',nph,nph0
c      print *,'entering depcor'
      if(nph.eq.nph0) go to 1
      nph0=nph
      us(nph)=umod(zs,isrc,nph)
c   If we are in a high slowness zone, find the slowness of the lid.
      umin=us(nph)
      ks=isrc(nph)
c     write(10,*)'ks us',ks,sngl(umin)
      do 2 i=1,ks
      if(pm(i,nph).gt.umin) go to 2
      umin=pm(i,nph)
 2    continue
c   Find where the source slowness falls in the ray parameter array.
      n1=ku(nph)+1
      do 3 i=2,n1
      if(pu(i,nph).gt.umin) go to 4
 3    continue
      k2=n1
      if(pu(n1,nph).eq.umin) go to 50
      write(6,*)'Source slowness too large.'
      call abort
 4    k2=i
c50   write(10,*)'k2 umin',k2,sngl(umin)
c
c   Read in the appropriate depth correction values.
c
 50   noext=.false.
      sgn=1d0
      if(msrc(nph).eq.0) msrc(nph)=1
c   See if the source depth coincides with a model sample
      ztol=xn*tol/(1.-xn*odep)
      if(dabs(zs-zm(ks+1,nph)).gt.ztol) go to 5
      ks=ks+1
      go to 6
 5    if(dabs(zs-zm(ks,nph)).gt.ztol) go to 7
c   If so flag the fact and make sure that the right integrals are
c   available.
 6    noext=.true.
      if(msrc(nph).eq.ks) go to 8
      call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
      go to 11
c   If it is necessary to interpolate, see if appropriate integrals
c   have already been read in.
 7    if(msrc(nph).ne.ks+1) go to 9
      ks=ks+1
      sgn=-1d0
      go to 8
 9    if(msrc(nph).eq.ks) go to 8
c   If not, read in integrals for the model depth nearest the source
c   depth.
      if(dabs(zm(ks,nph)-zs).le.dabs(zm(ks+1,nph)-zs)) go to 10
      ks=ks+1
      sgn=-1d0
 10   call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
c   Move the depth correction values to a less temporary area.
 11   do 31 i=1,ku(nph)
 31   tauu(i,nph)=tup(i)
      k=ku(nph)
      do 12 i=1,km(nph)
      k=k+1
      xc(i)=tup(k)
 12   xu(i,nph)=tup(k)
c     write(10,*)'bkin',ks,sngl(sgn),sngl(tauu(1,nph)),sngl(xu(1,nph))
c
c   Fiddle pointers.
c
 8    msrc(nph)=ks
c     write(10,*)'msrc sgn',msrc(nph),sngl(sgn)
      noend=.false.
      if(dabs(umin-pu(k2-1,nph)).le.dtol*umin) k2=k2-1
      if(dabs(umin-pu(k2,nph)).le.dtol*umin) noend=.true.
      if(msrc(nph).le.1.and.noext) msrc(nph)=0
      k1=k2-1
      if(noend) k1=k2
c     write(10,*)'noend noext k2 k1',noend,noext,k2,k1
      if(noext) go to 14
c
c   Correct the integrals for the depth interval [zm(msrc),zs].
c
      ms=msrc(nph)
      if(sgn)15,16,16
 16   u0=pm(ms,nph)
      z0=zm(ms,nph)
      u1=us(nph)
      z1=zs
      go to 17
 15   u0=us(nph)
      z0=zs
      u1=pm(ms,nph)
      z1=zm(ms,nph)
 17   mu=1
c     write(10,*)'u0 z0',sngl(u0),sngl(z0)
c     write(10,*)'u1 z1',sngl(u1),sngl(z1)
      do 18 k=1,k1
      call tauint(pu(k,nph),u0,u1,z0,z1,ttau,tx)
      tauc(k)=tauu(k,nph)+sgn*ttau
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 18
      xc(mu)=xu(mu,nph)+sgn*tx
c     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 18   continue
      go to 39
c   If there is no correction, copy the depth corrections to working
c   storage.
 14   mu=1
      do 40 k=1,k1
      tauc(k)=tauu(k,nph)
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 40
      xc(mu)=xu(mu,nph)
c     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 40   continue
c
c   Calculate integrals for the ray bottoming at the source depth.
c
 39   xus1(nph)=0d0
      xus2(nph)=0d0
      mu=mu-1
      if(dabs(umin-us(nph)).gt.dtol.and.dabs(umin-pux(mu,nph)).le.dtol)
     1  mu=mu-1
c   This loop may be skipped only for surface focus as range is not
c   available for all ray parameters.
      if(msrc(nph).le.0) go to 1
      is=isrc(nph)
      tauus2(nph)=0d0
      if(dabs(pux(mu,nph)-umin).gt.dtol.or.dabs(us(nph)-umin).gt.dtol)
     1  go to 48
c   If we happen to be right at a discontinuity, range is available.
      tauus1(nph)=tauc(k1)
      xus1(nph)=xc(mu)
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph)),'  *'
      go to 33
c   Integrate from the surface to the source.
 48   tauus1(nph)=0d0
      j=1
      if(is.lt.2) go to 42
      do 19 i=2,is
      call tauint(umin,pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
 19   j=i
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph))
 42   if(dabs(zm(is,nph)-zs).le.dtol) go to 33
c   Unless the source is right on a sample slowness, one more partial
c   integral is needed.
      call tauint(umin,pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph))
 33   if(pm(is+1,nph).lt.umin) go to 41
c   If we are in a high slowness zone, we will also need to integrate
c   down to the turning point of the shallowest down-going ray.
      u1=us(nph)
      z1=zs
      do 35 i=is+1,mt(nph)
      u0=u1
      z0=z1
      u1=pm(i,nph)
      z1=zm(i,nph)
      if(u1.lt.umin) go to 36
      call tauint(umin,u0,u1,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
 35   xus2(nph)=xus2(nph)+tx
c36   write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)),
c    1 sngl(xus2(nph))
 36   z1=zmod(umin,i-1,nph)
      if(dabs(z0-z1).le.dtol) go to 41
c   Unless the turning point is right on a sample slowness, one more
c   partial integral is needed.
      call tauint(umin,u0,umin,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
      xus2(nph)=xus2(nph)+tx
c     write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)),
c    1 sngl(xus2(nph))
c
c   Take care of converted phases.
c
 41   iph=mod(nph,2)+1
      xus1(iph)=0d0
      xus2(iph)=0d0
      tauus1(iph)=0d0
      tauus2(iph)=0d0
      go to (59,61),nph
 61   if(umin.gt.pu(ku(1)+1,1)) go to 53
c
c   If we are doing an S-wave depth correction, we may need range and
c   tau for the P-wave which turns at the S-wave source slowness.  This
c   would bd needed for sPg and SPg when the source is in the deep mantle.
c
      do 44 j=1,nbrn
      if((phcd(j)(1:2).ne.'sP'.and.phcd(j)(1:2).ne.'SP').or.
     1 px(j,2).le.0d0) go to 44
c     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1),
c    1 px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 44   continue
      go to 53
c
c   If we are doing an P-wave depth correction, we may need range and
c   tau for the S-wave which turns at the P-wave source slowness.  This
c   would be needed for pS and PS.
c
 59   do 60 j=1,nbrn
      if((phcd(j)(1:2).ne.'pS'.and.phcd(j)(1:2).ne.'PS').or.
     1 px(j,2).le.0d0) go to 60
c     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1),
c    1 px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 60   continue
      go to 53
c
c   Do the integral.
 45   j=1
c     write(10,*)'Depcor:  do pS or sP integral - iph =',iph
      do 46 i=2,mt(iph)
      if(umin.ge.pm(i,iph)) go to 47
      call tauint(umin,pm(j,iph),pm(i,iph),zm(j,iph),zm(i,iph),ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
 46   j=i
 47   z1=zmod(umin,j,iph)
      if(dabs(zm(j,iph)-z1).le.dtol) go to 53
c   Unless the turning point is right on a sample slowness, one more
c   partial integral is needed.
      call tauint(umin,pm(j,iph),umin,zm(j,iph),z1,ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
c     write(10,*)'is ks tauusp xusp',j,ks,sngl(tauus1(iph)),
c    1 sngl(xus1(iph))
c
 53   ua(1,nph)=-1d0
c     if(odep.ge.deplim.or.odep.le..1) go to 43
      if(odep.ge.deplim) go to 43
      do 57 i=1,nseg
      if(.not.segmsk(i)) go to 57
      if(nafl(i,1).eq.nph.and.nafl(i,2).eq.0.and.iidx(i).le.0) go to 58
 57   continue
      go to 43
c
c   If the source is very shallow, we will need to insert some extra
c   ray parameter samples into the up-going branches.
c
 58   du=amin1(1e-5+(odep-.4)*2e-5,1e-5)
c     write(10,*)'Add:  nph is ka odep du us =',nph,is,ka,odep,
c    1 sngl(du),sngl(us(nph))
      lp=lpower
      k=0
      do 56 l=ka,1,-1
      k=k+1
      ua(k,nph)=us(nph)-(l**lp)*du
      lp=lp-1
      taua(k,nph)=0d0
      j=1
      if(is.lt.2) go to 54
      do 55 i=2,is
      call tauint(ua(k,nph),pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),
     1 ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
 55   j=i
c     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 54   if(dabs(zm(is,nph)-zs).le.dtol) go to 56
c   Unless the source is right on a sample slowness, one more partial
c   integral is needed.
      call tauint(ua(k,nph),pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
c     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 56   continue
      go to 43
c
c   Construct tau for all branches.
c

 1    mu=mu+1
 43   j=1
c     write(10,*)'mu',mu
c     write(10,*)'killer loop:'
      do 20 i=1,nseg
      if(.not.segmsk(i)) go to 20
c     write(10,*)'i iidx nafl nph',i,iidx(i),nafl(i,1),nph
      if(iidx(i).gt.0.or.iabs(nafl(i,1)).ne.nph.or.(msrc(nph).le.0.and.
     1 nafl(i,1).gt.0)) go to 20
c
      iph=nafl(i,2)
      kph=nafl(i,3)
c   Handle up-going P and S.
      if(iph.le.0) iph=nph
      if(kph.le.0) kph=nph
      sgn=isign(1,nafl(i,1))
      i1=indx(i,1)
      i2=indx(i,2)
c     write(10,*)'i1 i2 sgn iph',i1,i2,sngl(sgn),iph
      m=1
      do 21 k=i1,i2
      if(pt(k).gt.umin) go to 22
 23   if(dabs(pt(k)-pu(m,nph)).le.dtol) go to 2115
      m=m+1
      go to 23
 2115 tau(1,k)=taut(k)+sgn*tauc(m)
 21   continue 

      k=i2
c     write(10,*)'k m',k,m
      go to 24
c22   write(10,*)'k m',k,m
 22   if(dabs(pt(k-1)-umin).le.dtol) k=k-1
      ki=ki+1
      kk(ki)=k
      pk(ki)=pt(k)
      pt(k)=umin
      fac=fcs(i,1)
c     write(10,*)'ki fac',ki,sngl(fac)
      tau(1,k)=fac*(tauus1(iph)+tauus2(iph)+tauus1(kph)+tauus2(kph))+
     1 sgn*tauus1(nph)
c     write(10,*)'&&&&& nph iph kph tauus1 tauus2 tau =',
c    1 nph,iph,kph,sngl(tauus1(1)),sngl(tauus1(2)),sngl(tauus2(1)),
c    2 sngl(tauus2(2)),sngl(tau(1,k))
 24   m=1
 26   if(jndx(j,1).ge.indx(i,1)) go to 25
      j=j+1
      go to 26
 25   jndx(j,2)=min0(jidx(j),k)
      if(jndx(j,1).lt.jndx(j,2)) go to 37
      jndx(j,2)=-1
      go to 20
c37   write(10,*)'j jndx jidx',j,jndx(j,1),jndx(j,2),jidx(j),' ',
c    1 phcd(j)
 37   do 30 l=1,2
 28   if(dabs(pux(m,nph)-px(j,l)).le.dtol) go to 27
      if(m.ge.mu) go to 29
      m=m+1
      go to 28
 27   xbrn(j,l)=xt(j,l)+sgn*xc(m)
c     write(10,*)'x up:  j l m  ',j,l,m
      go to 30
 29   xbrn(j,l)=fac*(xus1(iph)+xus2(iph)+xus1(kph)+xus2(kph))+
     1 sgn*xus1(nph)
c     write(10,*)'x up:  j l end',j,l
c     write(10,*)'&&&&& nph iph kph xus1 xus2 xbrn =',
c    1 nph,iph,kph,sngl(xus1(1)),sngl(xus1(2)),sngl(xus2(1)),
c    2 sngl(xus2(2)),sngl(xbrn(j,l))
 30   continue
      if(j.ge.nbrn) go to 20
      j=j+1
      if(jndx(j,1).le.k) go to 25
 20   continue
      return
      end

      double precision function umod(zs,isrc,nph)
      save 
      include 'ttlim.inc'
      character*31 msg
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn
      double precision zs,uend,dtol,zmod
      dimension isrc(2)
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      data dtol/1d-6/
c
      m1=mt(nph)
      do 1 i=2,m1
      if(zm(i,nph).le.zs) go to 2
 1    continue
      dep=(1d0-dexp(zs))/xn
      write(msg,100)dep
      write(6,100)dep
 100  format('Source depth (',f6.1,') too deep.')
      write(6,*)msg
      call abort
 2    if(dabs(zs-zm(i,nph)).le.dtol.and.dabs(zm(i,nph)-zm(i+1,nph)).le.
     1 dtol) go to 3
      j=i-1
      isrc(nph)=j
      umod=pm(j,nph)+(pm(i,nph)-pm(j,nph))*(dexp(zs-zm(j,nph))-1d0)/
     1 (dexp(zm(i,nph)-zm(j,nph))-1d0)
      return
 3    isrc(nph)=i
      umod=pm(i+1,nph)
      return
c
      entry zmod(uend,js,nph)
      i=js+1
      zmod=zm(js,nph)+dlog(dmax1((uend-pm(js,nph))*(dexp(zm(i,nph)-
     1 zm(js,nph))-1d0)/(pm(i,nph)-pm(js,nph))+1d0,1d-30))
      return
      end

      subroutine bkin(lu,nrec,len,buf)
c
c $$$$$ calls no other routines $$$$$
c
c   Bkin reads a block of len double precision words into array buf(len)
c   from record nrec of the direct access unformatted file connected to
c   logical unit lu.
c
      save
      double precision buf(len),tmp
c
      if(nrec.le.0) go to 1
      read(lu,rec=nrec)buf
      tmp=buf(1)
      return
c   If the record doesn't exist, zero fill the buffer.
 1    do 2 i=1,len
 2    buf(i)=0d0
      return
      end

      subroutine tauint(ptk,ptj,pti,zj,zi,tau,x)
      save
c
c $$$$$ calls warn $$$$$
c
c   Tauint evaluates the intercept (tau) and distance (x) integrals  for
c   the spherical earth assuming that slowness is linear between radii
c   for which the model is known.  The partial integrals are performed
c   for ray slowness ptk between model radii with slownesses ptj and pti
c   with equivalent flat earth depths zj and zi respectively.  The partial
c   integrals are returned in tau and x.  Note that ptk, ptj, pti, zj, zi,
c   tau, and x are all double precision.
c
      character*71 msg
      double precision ptk,ptj,pti,zj,zi,tau,x
      double precision xx,b,sqk,sqi,sqj,sqb
c
      if(dabs(zj-zi).le.1d-9) go to 13
      if(dabs(ptj-pti).gt.1d-9) go to 10
      if(dabs(ptk-pti).le.1d-9) go to 13
      b=dabs(zj-zi)
      sqj=dsqrt(dabs(ptj*ptj-ptk*ptk))
      tau=b*sqj
      x=b*ptk/sqj
      go to 4
 10   if(ptk.gt.1d-9.or.pti.gt.1d-9) go to 1
c   Handle the straight through ray.
      tau=ptj
      x=1.5707963267948966d0
      go to 4
 1    b=ptj-(pti-ptj)/(dexp(zi-zj)-1d0)
      if(ptk.gt.1d-9) go to 2
      tau=-(pti-ptj+b*dlog(pti/ptj)-b*dlog(dmax1((ptj-b)*pti/
     1 ((pti-b)*ptj),1d-30)))
      x=0d0
      go to 4
 2    if(ptk.eq.pti) go to 3
      if(ptk.eq.ptj) go to 11
      sqk=ptk*ptk
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(sqb.gt.1d-30) go to 5
      xx=0d0
      x=ptk*(dsqrt(dabs((pti+b)/(pti-b)))-dsqrt(dabs((ptj+b)/
     1 (ptj-b))))/b
      go to 6
 5    if(b*b.lt.sqk) go to 7
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*
     1 (sqb*sqj+b*ptj-sqk)),1d-30))
      x=ptk*xx/sqb
      go to 6
 7    xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptk*dabs(pti-b)),1d0),-1d0))-
     1 dasin(dmax1(dmin1((b*ptj-sqk)/(ptk*dabs(ptj-b)),1d0),-1d0))
      x=-ptk*xx/sqb
 6    tau=-(sqi-sqj+b*dlog((pti+sqi)/(ptj+sqj))-sqb*xx)
      go to 4
 3    sqk=pti*pti
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 8
      xx=dlog(dmax1((ptj-b)*(b*pti-sqk)/((pti-b)*(sqb*sqj+b*ptj-sqk)),
     1 1d-30))
      x=pti*xx/sqb
      go to 9
 8    xx=dsign(1.5707963267948966d0,b-pti)-dasin(dmax1(dmin1((b*ptj-
     1 sqk)/(pti*dabs(ptj-b)),1d0),-1d0))
      x=-pti*xx/sqb
 9    tau=-(b*dlog(pti/(ptj+sqj))-sqj-sqb*xx)
      go to 4
 11   sqk=ptj*ptj
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 12
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*(b*ptj-sqk)),
     1 1d-30))
      x=ptj*xx/sqb
      go to 14
 12   xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptj*dabs(pti-b)),1d0),-1d0))-
     1 dsign(1.5707963267948966d0,b-ptj)
      x=-ptj*xx/sqb
 14   tau=-(b*dlog((pti+sqi)/ptj)+sqi-sqb*xx)
c
c   Handle various error conditions.
c
 4    if(x.ge.-1d-10) go to 15
      write(msg,100)ptk,ptj,pti,tau,x
 100  format('Bad range: ',1p5d12.4)
      call warn(msg)
 15   if(tau.ge.-1d-10) go to 16
      write(msg,101)ptk,ptj,pti,tau,x
 101  format('Bad tau: ',1p5d12.4)
      call warn(msg(1:69))
 16   return
c   Trap null integrals and handle them properly.
 13   tau=0d0
      x=0d0
      return
      end


      subroutine spfit(jb,int)
      save 
      include 'ttlim.inc'
      character*3 disc
      character*8 phcd
      logical newgrd,makgrd,segmsk,prnt
c     logical log
      double precision pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
      double precision pmn,dmn,dmx,hm,shm,thm,p0,p1,tau0,tau1,x0,x1,pe,
     1 pe0,spe0,scpe0,pe1,spe1,scpe1,dpe,dtau,dbrnch,cn,x180,x360,dtol,
     2 ptol,xmin
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/prtflc/segmsk(jseg),prnt(2)
      data dbrnch,cn,x180,x360,xmin,dtol,ptol/2.5307274d0,57.295779d0,
     1 3.1415927d0,6.2831853d0,3.92403d-3,1d-6,2d-6/
c
      if(prnt(1)) write(10,102)
      i1=jndx(jb,1)
      i2=jndx(jb,2)
c     write(10,*)'Spfit:  jb i1 i2 pt =',jb,i1,i2,sngl(pt(i1)),
c    1 sngl(pt(i2))
      if(i2-i1.gt.1.or.dabs(pt(i2)-pt(i1)).gt.ptol) go to 14
      jndx(jb,2)=-1
      return
 14   newgrd=.false.
      makgrd=.false.
      if(dabs(px(jb,2)-pt(i2)).gt.dtol) newgrd=.true.
c     write(10,*)'Spfit:  px newgrd =',sngl(px(jb,2)),newgrd
      if(.not.newgrd) go to 10
      k=mod(int-1,2)+1
      if(int.ne.int0(k)) makgrd=.true.
c     write(10,*)'Spfit:  int k int0 makgrd =',int,k,int0(k),makgrd
      if(int.gt.2) go to 12
c     call query('Enter xmin:',log)
c     read *,xmin
c     xmin=xmin*xn
      xmin=xn*amin1(amax1(2.*odep,2.),25.)
c     write(10,*)'Spfit:  xmin =',xmin,xmin/xn
      call pdecu(i1,i2,xbrn(jb,1),xbrn(jb,2),xmin,int,i2)
      jndx(jb,2)=i2
 12   nn=i2-i1+1
      if(makgrd) call tauspl(1,nn,pt(i1),tcoef(1,1,k))
c     write(10,301,iostat=ios)jb,k,nn,int,newgrd,makgrd,
c    1 xbrn(jb,1),xbrn(jb,2),(i,pt(i-1+i1),tau(1,i-1+i1),
c    2 (tcoef(j,i,k),j=1,5),i=1,nn)
c301  format(/1x,4i3,2l3,2f12.8/(1x,i5,0p2f12.8,1p5d10.2))
      call fitspl(1,nn,tau(1,i1),xbrn(jb,1),xbrn(jb,2),tcoef(1,1,k))
      int0(k)=int
      go to 11
 10   call fitspl(i1,i2,tau,xbrn(jb,1),xbrn(jb,2),coef)
 11   pmn=pt(i1)
      dmn=xbrn(jb,1)
      dmx=dmn
      mxcnt=0
      mncnt=0
c     call appx(i1,i2,xbrn(jb,1),xbrn(jb,2))
c     write(10,300)(i,pt(i),(tau(j,i),j=1,3),i=i1,i2)
c300  format(/(1x,i5,4f12.6))
      pe=pt(i2)
      p1=pt(i1)
      tau1=tau(1,i1)
      x1=tau(2,i1)
      pe1=pe-p1
      spe1=dsqrt(dabs(pe1))
      scpe1=pe1*spe1
      j=i1
      is=i1+1
      do 2 i=is,i2
      p0=p1
      p1=pt(i)
      tau0=tau1
      tau1=tau(1,i)
      x0=x1
      x1=tau(2,i)
      dpe=p0-p1
      dtau=tau1-tau0
      pe0=pe1
      pe1=pe-p1
      spe0=spe1
      spe1=dsqrt(dabs(pe1))
      scpe0=scpe1
      scpe1=pe1*spe1
      tau(4,j)=(2d0*dtau-dpe*(x1+x0))/(.5d0*(scpe1-scpe0)-1.5d0*spe1*
     1 spe0*(spe1-spe0))
      tau(3,j)=(dtau-dpe*x0-(scpe1+.5d0*scpe0-1.5d0*pe1*spe0)*tau(4,j))/
     1 (dpe*dpe)
      tau(2,j)=(dtau-(pe1*pe1-pe0*pe0)*tau(3,j)-(scpe1-scpe0)*tau(4,j))/
     1 dpe
      tau(1,j)=tau0-scpe0*tau(4,j)-pe0*(pe0*tau(3,j)+tau(2,j))
      xlim(1,j)=dmin1(x0,x1)
      xlim(2,j)=dmax1(x0,x1)
      if(xlim(1,j).ge.dmn) go to 5
      dmn=xlim(1,j)
      pmn=pt(j)
      if(x1.lt.x0) pmn=pt(i)
 5    disc=' '
      if(dabs(tau(3,j)).le.1d-30) go to 4
      shm=-.375d0*tau(4,j)/tau(3,j)
      hm=shm*shm
      if(shm.le.0d0.or.(hm.le.pe1.or.hm.ge.pe0)) go to 4
 7    thm=tau(2,j)+shm*(2d0*shm*tau(3,j)+1.5d0*tau(4,j))
      xlim(1,j)=dmin1(xlim(1,j),thm)
      xlim(2,j)=dmax1(xlim(2,j),thm)
      if(thm.ge.dmn) go to 6
      dmn=thm
      pmn=pe-hm
 6    disc='max'
      if(tau(4,j).lt.0d0) disc='min'
      if(disc.eq.'max') mxcnt=mxcnt+1
      if(disc.eq.'min') mncnt=mncnt+1
 4    if(prnt(1)) write(10,100,iostat=ios)disc,j,pt(j),
     1 (tau(k,j),k=1,4),(cn*xlim(k,j),k=1,2)
 100  format(1x,a,i5,f10.6,1p4e10.2,0p2f7.2)
      dmx=dmax1(dmx,xlim(2,j))
 2    j=i
c     if(prnt(1)) write(10,100,iostat=ios)'   ',j,pt(j)
      xbrn(jb,1)=dmn
      xbrn(jb,2)=dmx
      xbrn(jb,3)=pmn
      idel(jb,1)=1
      idel(jb,2)=1
      if(xbrn(jb,1).gt.x180) idel(jb,1)=2
      if(xbrn(jb,2).gt.x180) idel(jb,2)=2
      if(xbrn(jb,1).gt.x360) idel(jb,1)=3
      if(xbrn(jb,2).gt.x360) idel(jb,2)=3
      if(int.gt.2) go to 1
      phcd(jb)=phcd(jb)(1:1)
      i=jb
      do 8 j=1,nbrn
      i=mod(i,nbrn)+1
      if(phcd(i)(1:1).eq.phcd(jb).and.phcd(i)(2:2).ne.'P'.and.
     1 (pe.ge.px(i,1).and.pe.le.px(i,2))) go to 9
 8    continue
      go to 1
 9    phcd(jb)=phcd(i)
      if(dabs(pt(i2)-pt(jndx(i,1))).le.dtol) phcd(jb)=phcd(i-1)
 1    if(prnt(1).and.prnt(2)) write(10,102)
 102  format()
      if(dbrn(jb,1).le.0d0) go to 3
      dbrn(jb,1)=dmx
      dbrn(jb,2)=dbrnch
      if(prnt(2)) write(10,101,iostat=ios)phcd(jb),
     1 (jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),
     2 (cn*dbrn(jb,k),k=1,2),(idel(jb,k),k=1,3),int,newgrd,makgrd
 101  format(1x,a,2i5,2f8.2,f8.4,2f8.2,4i3,2l2)
      go to 15
 3    if(prnt(2)) write(10,103,iostat=ios)phcd(jb),
     1 (jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),
     2 (idel(jb,k),k=1,3),int,newgrd,makgrd
 103  format(1x,a,2i5,2f8.2,f8.4,16x,4i3,2l2)
 15   if(mxcnt.gt.mncnt.or.mncnt.gt.mxcnt+1)
     1 call warn('Bad interpolation on '//phcd(jb))
      return
      end
c
c-----------------------------------------------------------------
c
      subroutine pdecu(i1,i2,x0,x1,xmin,int,len)
      save 
      include 'ttlim.inc'
      double precision us,pt,tau,xlim,xbrn,dbrn,ua,taua
      double precision x0,x1,xmin,dx,dx2,sgn,rnd,xm,axm,x,h1,h2,hh,xs
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/pdec/ua(5,2),taua(5,2),deplim,ka
c
c     write(10,*)'Pdecu:  ua =',sngl(ua(1,int))
      if(ua(1,int).le.0d0) go to 17
c     write(10,*)'Pdecu:  fill in new grid'
      k=i1+1
      do 18 i=1,ka
      pt(k)=ua(i,int)
      tau(1,k)=taua(i,int)
 18   k=k+1
      pt(k)=pt(i2)
      tau(1,k)=tau(1,i2)
      go to 19
c
 17   is=i1+1
      ie=i2-1
      xs=x1
      do 11 i=ie,i1,-1
      x=xs
      if(i.ne.i1) go to 12
      xs=x0
      go to 14
 12   h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      xs=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 14   if(dabs(x-xs).le.xmin) go to 15
 11   continue
      len=i2
      return
 15   ie=i
      if(dabs(x-xs).gt..75d0*xmin.or.ie.eq.i2) go to 16
      xs=x
      ie=ie+1
 16   n=max0(idint(dabs(xs-x0)/xmin+.8d0),1)
      dx=(xs-x0)/n
      dx2=dabs(.5d0*dx)
      sgn=dsign(1d0,dx)
      rnd=0d0
      if(sgn.gt.0d0) rnd=1d0
      xm=x0+dx
      k=i1
      m=is
      axm=1d10
      do 1 i=is,ie
      if(i.lt.ie) go to 8
      x=xs
      go to 5
 8    h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      x=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 5    if(sgn*(x-xm).le.dx2) go to 2
      if(k.lt.m) go to 3
      do 4 j=m,k
 4    pt(j)=-1d0
 3    m=k+2
      k=i-1
      axm=1d10
 7    xm=xm+dx*idint((x-xm-dx2)/dx+rnd)
 2    if(dabs(x-xm).ge.axm) go to 1
      axm=dabs(x-xm)
      k=i-1
 1    continue
      if(k.lt.m) go to 9
      do 6 j=m,k
 6    pt(j)=-1d0
 9    k=i1
      do 10 i=is,i2
      if(pt(i).lt.0d0) go to 10
      k=k+1
      pt(k)=pt(i)
      tau(1,k)=tau(1,i)
 10   continue
 19   len=k
c     write(10,300)(i,pt(i),tau(1,i),i=i1,len)
c300  format(/(1x,i5,0pf12.6,1pd15.4))
      return
      end
      subroutine tauspl(i1,i2,pt,coef)
c
c $$$$$ calls only library routines $$$$$
c
c   Given ray parameter grid pt;i (pt sub i), i=i1,i1+1,...,i2, tauspl
c   determines the i2-i1+3 basis functions for interpolation I such
c   that:
c
c      tau(p) = a;1,i + Dp * a;2,i + Dp**2 * a;3,i + Dp**(3/2) * a;4,i
c
c   where Dp = pt;n - p, pt;i <= p < pt;i+1, and the a;j,i's are
c   interpolation coefficients.  Rather than returning the coefficients,
c   a;j,i, which necessarily depend on tau(pt;i), i=i1,i1+1,...,i2 and
c   x(pt;i) (= -d tau(p)/d p | pt;i), i=i1,i2, tauspl returns the
c   contribution of each basis function and its derivitive at each
c   sample.  Each basis function is non-zero at three grid points,
c   therefore, each grid point will have contributions (function values
c   and derivitives) from three basis functions.  Due to the basis
c   function normalization, one of the function values will always be
c   one and is not returned in array coef with the other values.
c   Rewritten on 23 December 1983 by R. Buland.
c
      save
      double precision pt(i2),coef(5,i2)
      double precision del(5),sdel(5),deli(5),d3h(4),d1h(4),dih(4),
     1 d(4),ali,alr,b3h,b1h,bih,th0p,th2p,th3p,th2m
c
      n2=i2-i1-1
      if(n2.le.-1) return
      is=i1+1
c
c   To achieve the requisite stability, proceed by constructing basis
c   functions G;i, i=0,1,...,n+1.  G;i will be non-zero only on the
c   interval [p;i-2,p;i+2] and will be continuous with continuous first
c   and second derivitives.  G;i(p;i-2) and G;i(p;i+2) are constrained
c   to be zero with zero first and second derivitives.  G;i(p;i) is
c   normalized to unity.
c
c   Set up temporary variables appropriate for G;-1.  Note that to get
c   started, the ray parameter grid is extrapolated to yeild p;i, i=-2,
c   -1,0,1,...,n.
      del(2)=pt(i2)-pt(i1)+3d0*(pt(is)-pt(i1))
      sdel(2)=dsqrt(dabs(del(2)))
      deli(2)=1d0/sdel(2)
      m=2
      do 1 k=3,5
      del(k)=pt(i2)-pt(i1)+(5-k)*(pt(is)-pt(i1))
      sdel(k)=dsqrt(dabs(del(k)))
      deli(k)=1d0/sdel(k)
      d3h(m)=del(k)*sdel(k)-del(m)*sdel(m)
      d1h(m)=sdel(k)-sdel(m)
      dih(m)=deli(k)-deli(m)
 1    m=k
      l=i1-1
      if(n2.le.0) go to 10
c   Loop over G;i, i=0,1,...,n-3.
      do 2 i=1,n2
      m=1
c   Update temporary variables for G;i-1.
      do 3 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 3
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 3    m=k
      l=l+1
      del(5)=pt(i2)-pt(l+1)
      sdel(5)=dsqrt(dabs(del(5)))
      deli(5)=1d0/sdel(5)
      d3h(4)=del(5)*sdel(5)-del(4)*sdel(4)
      d1h(4)=sdel(5)-sdel(4)
      dih(4)=deli(5)-deli(4)
c   Construct G;i-1.
      ali=1d0/(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*dih(1)*del(3))*
     1 del(3))
      alr=ali*(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*del(3)*
     1 deli(2)-sdel(3))*del(3))
      b3h=d3h(2)+alr*d3h(1)
      b1h=d1h(2)+alr*d1h(1)
      bih=dih(2)+alr*dih(1)
      th0p=d1h(1)*b3h-d3h(1)*b1h
      th2p=d1h(3)*b3h-d3h(3)*b1h
      th3p=d1h(4)*b3h-d3h(4)*b1h
      th2m=dih(3)*b3h-d3h(3)*bih
c   The d;i's completely define G;i-1.
      d(4)=ali*((dih(1)*b3h-d3h(1)*bih)*th2p-th2m*th0p)/((dih(4)*b3h-
     1 d3h(4)*bih)*th2p-th2m*th3p)
      d(3)=(th0p*ali-th3p*d(4))/th2p
      d(2)=(d3h(1)*ali-d3h(3)*d(3)-d3h(4)*d(4))/b3h
      d(1)=alr*d(2)-ali
c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
c   G;i-1(p;i-1) need not be constructed as it is normalized to unity.
      coef(1,l)=(.125d0*del(5)*sdel(5)-(.75d0*sdel(5)+.375d0*deli(5)*
     1 del(4)-sdel(4))*del(4))*d(4)
      if(i.ge.3) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+
     1 .375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      coef(3,l)=-.75d0*(sdel(5)+deli(5)*del(4)-2d0*sdel(4))*d(4)
      if(i.ge.2) coef(4,l-1)=-.75d0*((sdel(2)+deli(2)*del(3)-
     1 2d0*sdel(3))*d(2)-(d1h(1)+dih(1)*del(3))*d(1))
      if(i.ge.3) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-
     1 2d0*sdel(2))*d(1)
 2    continue
c   Loop over G;i, i=n-2,n-1,n,n+1.  These cases must be handled
c   seperately because of the singularities in the second derivitive
c   at p;n.
 10   do 4 j=1,4
      m=1
c   Update temporary variables for G;i-1.
      do 5 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 5
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 5    m=k
      l=l+1
      del(5)=0d0
      sdel(5)=0d0
      deli(5)=0d0
c   Construction of the d;i's is different for each case.  In cases
c   G;i, i=n-1,n,n+1, G;i is truncated at p;n to avoid patching across
c   the singularity in the second derivitive.
      if(j.lt.4) go to 6
c   For G;n+1 constrain G;n+1(p;n) to be .25.
      d(1)=2d0/(del(1)*sdel(1))
      go to 9
c   For G;i, i=n-2,n-1,n, the condition dG;i(p)/dp|p;i = 0 has been
c   substituted for the second derivitive continuity condition that
c   can no longer be satisfied.
 6    alr=(sdel(2)+deli(2)*del(3)-2d0*sdel(3))/(d1h(1)+dih(1)*del(3))
      d(2)=1d0/(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*deli(2)*
     1 del(3)-sdel(3))*del(3)-(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*
     2 dih(1)*del(3))*del(3))*alr)
      d(1)=alr*d(2)
      if(j-2)8,7,9
c   For G;n-1 constrain G;n-1(p;n) to be .25.
 7    d(3)=(2d0+d3h(2)*d(2)+d3h(1)*d(1))/(del(3)*sdel(3))
      go to 9
c   No additional constraints are required for G;n-2.
 8    d(3)=-((d3h(2)-d1h(2)*del(4))*d(2)+(d3h(1)-d1h(1)*del(4))*
     1 d(1))/(d3h(3)-d1h(3)*del(4))
      d(4)=(d3h(3)*d(3)+d3h(2)*d(2)+d3h(1)*d(1))/(del(4)*sdel(4))
c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
 9    if(j.le.2) coef(1,l)=(.125d0*del(3)*sdel(3)-(.75d0*sdel(3)+.375d0*
     1 deli(3)*del(4)-sdel(4))*del(4))*d(3)-(.125d0*d3h(2)-(.75d0*
     2 d1h(2)+.375d0*dih(2)*del(4))*del(4))*d(2)-(.125d0*d3h(1)-(.75d0*
     3 d1h(1)+.375d0*dih(1)*del(4))*del(4))*d(1)
      if(l-i1.gt.1) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+
     1 .375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      if(j.le.2) coef(3,l)=-.75d0*((sdel(3)+deli(3)*del(4)-
     1 2d0*sdel(4))*d(3)-(d1h(2)+dih(2)*del(4))*d(2)-(d1h(1)+
     2 dih(1)*del(4))*d(1))
      if(j.le.3.and.l-i1.gt.0) coef(4,l-1)=0d0
      if(l-i1.gt.1) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-
     1 2d0*sdel(2))*d(1)
 4    continue
      return
      end
      subroutine fitspl(i1,i2,tau,x1,xn,coef)
c
c $$$$$ calls only library routines $$$$$
c
c   Given ray parameter grid p;i (p sub i), i=1,2,...,n, corresponding
c   tau;i values, and x;1 and x;n (x;i = -dtau/dp|p;i); tauspl finds
c   interpolation I such that:  tau(p) = a;1,i + Dp * a;2,i + Dp**2 *
c   a;3,i + Dp**(3/2) * a;4,i where Dp = p;n - p and p;i <= p < p;i+1.
c   Interpolation I has the following properties:  1) x;1, x;n, and
c   tau;i, i=1,2,...,n are fit exactly, 2) the first and second
c   derivitives with respect to p are continuous everywhere, and
c   3) because of the paramaterization d**2 tau/dp**2|p;n is infinite.
c   Thus, interpolation I models the asymptotic behavior of tau(p)
c   when tau(p;n) is a branch end due to a discontinuity in the
c   velocity model.  Note that array a must be dimensioned at least
c   a(4,n) though the interpolation coefficients will be returned in
c   the first n-1 columns.  The remaining column is used as scratch
c   space and returned as all zeros.  Programmed on 16 August 1982 by
c   R. Buland.
c
      save 
      double precision tau(4,i2),x1,xn,coef(5,i2),a(2,100),ap(3),
     1 b(100),alr,g1,gn
c
      if(i2-i1)13,1,2
 1    tau(2,i1)=x1
 13   return
 2    n=0
      do 3 i=i1,i2
      n=n+1
      b(n)=tau(1,i)
      do 3 j=1,2
 3    a(j,n)=coef(j,i)
      do 4 j=1,3
 4    ap(j)=coef(j+2,i2)
      n1=n-1
c
c   Arrays ap(*,1), a, and ap(*,2) comprise n+2 x n+2 penta-diagonal
c   matrix A.  Let x1, tau, and xn comprise corresponding n+2 vector b.
c   Then, A * g = b, may be solved for n+2 vector g such that
c   interpolation I is given by I(p) = sum(i=0,n+1) g;i * G;i(p).
c
c   Eliminate the lower triangular portion of A to form A'.  A
c   corresponding transformation applied to vector b is stored in
c   a(4,*).
      alr=a(1,1)/coef(3,i1)
      a(1,1)=1d0-coef(4,i1)*alr
      a(2,1)=a(2,1)-coef(5,i1)*alr
      b(1)=b(1)-x1*alr
      j=1
      do 5 i=2,n
      alr=a(1,i)/a(1,j)
      a(1,i)=1d0-a(2,j)*alr
      b(i)=b(i)-b(j)*alr
 5    j=i
      alr=ap(1)/a(1,n1)
      ap(2)=ap(2)-a(2,n1)*alr
      gn=xn-b(n1)*alr
      alr=ap(2)/a(1,n)
c   Back solve the upper triangular portion of A' for coefficients g;i.
c   When finished, storage g(2), a(4,*), g(5) will comprise vector g.
      gn=(gn-b(n)*alr)/(ap(3)-a(2,n)*alr)
      b(n)=(b(n)-gn*a(2,n))/a(1,n)
      j=n
      do 6 i=n1,1,-1
      b(i)=(b(i)-b(j)*a(2,i))/a(1,i)
 6    j=i
      g1=(x1-coef(4,i1)*b(1)-coef(5,i1)*b(2))/coef(3,i1)
c
      tau(2,i1)=x1
      is=i1+1
      ie=i2-1
      j=1
      do 7 i=is,ie
      j=j+1
 7    tau(2,i)=coef(3,i)*b(j-1)+coef(4,i)*b(j)+coef(5,i)*b(j+1)
      tau(2,i2)=xn
      return
      end
      subroutine trtm(delta,max,n,tt,dtdd,dtdh,dddp,phnm)
      save 
      include 'ttlim.inc'
      character*(*) phnm(max)
      character*8 ctmp(60)
      dimension tt(max),dtdd(max),dtdh(max),dddp(max),tmp(60,4),
     1 iptr(60)
      double precision us,pt,tau,xlim,xbrn,dbrn
      double precision x(3),cn,dtol,pi,pi2
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      data cn,dtol,atol,pi,pi2/.017453292519943296d0,1d-6,.005,
     1 3.1415926535897932d0,6.2831853071795865d0/
c
c      if (delta < 1.0) then
c      do i=1,jout
c         write(32,*) i,tau(1,i)
c      end do
c      stop
c      endif
      n=0
      if(mbr2.le.0) return
      x(1)=dmod(dabs(cn*delta),pi2)
      if(x(1).gt.pi) x(1)=pi2-x(1)
      x(2)=pi2-x(1)
      x(3)=x(1)+pi2
      if(dabs(x(1)).gt.dtol) go to 9
      x(1)=dtol
      x(3)=-10d0
 9    if(dabs(x(1)-pi).gt.dtol) go to 7
      x(1)=pi-dtol
      x(2)=-10d0
 7    do 1 j=mbr1,mbr2
 1    if(jndx(j,2).gt.0) call findtt(j,x,max,n,tmp,tmp(1,2),tmp(1,3),
     1 tmp(1,4),ctmp)
      if(n-1)3,4,5
 4    iptr(1)=1
      go to 6
 5    call r4sort(n,tmp,iptr)
 6    k=0
      do 2 i=1,n
      j=iptr(i)
      if(k.le.0) go to 8
      if(phnm(k).eq.ctmp(j).and.abs(tt(k)-tmp(j,1)).le.atol) go to 2
 8    k=k+1
      tt(k)=tmp(j,1)
      dtdd(k)=tmp(j,2)
      dtdh(k)=tmp(j,3)
      dddp(k)=tmp(j,4)
      phnm(k)=ctmp(j)
 2    continue
      n=k
 3    return
      end
      subroutine findtt(jb,x0,max,n,tt,dtdd,dtdh,dddp,phnm)
      save 
      include 'ttlim.inc'
      character*(*) phnm(max)
      character*8 phcd
      character*67 msg
      dimension tt(max),dtdd(max),dtdh(max),dddp(max)
      double precision us,pt,tau,xlim,xbrn,dbrn
      double precision x,x0(3),p0,p1,arg,dp,dps,delp,tol,ps,deps
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/pcdc/phcd(jbrn)
      data tol/3d-6/,deps/1d-10/
c
      nph=iabs(idel(jb,3))
      hsgn=isign(1,idel(jb,3))*hn
      dsgn=(-1.)**idel(jb,1)*dn
      dpn=-1./tn
      do 10 ij=idel(jb,1),idel(jb,2)
      x=x0(ij)
      dsgn=-dsgn
      if(x.lt.xbrn(jb,1).or.x.gt.xbrn(jb,2)) go to 12
      j=jndx(jb,1)
      is=j+1
      ie=jndx(jb,2)
      do 1 i=is,ie
      if(x.le.xlim(1,j).or.x.gt.xlim(2,j)) go to 8
      le=n
      p0=pt(ie)-pt(j)
      p1=pt(ie)-pt(i)
      delp=dmax1(tol*(pt(i)-pt(j)),1d-3)
      if(dabs(tau(3,j)).gt.1d-30) go to 2
      dps=(x-tau(2,j))/(1.5d0*tau(4,j))
      dp=dsign(dps*dps,dps)
      dp0=dp
      if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 9
      if(n.ge.max) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*(tau(1,j)+dp*(tau(2,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=dsgn*ps
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*.75d0*tau(4,j)/dmax1(dabs(dps),deps)
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 8
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
      go to 8
 2    do 4 jj=1,2
      go to (5,6),jj
 5    arg=9d0*tau(4,j)*tau(4,j)+32d0*tau(3,j)*(x-tau(2,j))
      if(arg.ge.0d0) go to 3
      write(msg,100)arg
 100  format('Bad sqrt argument:',1pd11.2,'.')
      call warn(msg(1:30))
 3    dps=-(3d0*tau(4,j)+dsign(dsqrt(dabs(arg)),tau(4,j)))/(8d0*
     1 tau(3,j))
      dp=dsign(dps*dps,dps)
      dp0=dp
      go to 7
 6    dps=(tau(2,j)-x)/(2d0*tau(3,j)*dps)
      dp=dsign(dps*dps,dps)
 7    if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 4
      if(n.ge.max) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=dsgn*ps
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*(2d0*tau(3,j)+.75d0*tau(4,j)/dmax1(dabs(dps),deps))
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 4
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
 4    continue
 9    if(n.gt.le) go to 8
      write(msg,101)phcd(jb),x,dp0,dp,p1,p0
 101  format('Failed to find phase:  ',a,f8.1,4f7.4)
      call warn(msg)
 8    j=i
 1    continue
c
 12   if(x.lt.dbrn(jb,1).or.x.gt.dbrn(jb,2)) go to 10
      if(n.ge.max) go to 13
      j=jndx(jb,1)
      i=jndx(jb,2)
      dp=pt(i)-pt(j)
      dps=dsqrt(dabs(dp))
      n=n+1
      tt(n)=tn*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+
     1 pt(j)*x)
      dtdd(n)=dsgn*sngl(pt(j))
      dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-pt(j)*pt(j))))
      dddp(n)=dpn*(2d0*tau(3,j)+.75d0*tau(4,j)/dmax1(dps,deps))
      ln=index(phcd(jb),' ')-1
      if(ln.le.0) ln=len(phcd(jb))
      phnm(n)=phcd(jb)(1:ln)//'diff'
 10   continue
      return
 13   write(msg,102)max
 102  format('More than',i3,' arrivals found.')
      call warn(msg(1:28))
      return
      end
      subroutine asnag1(lu,mode,n,ia,ib)
c
c $$$$$ calls assign, iargc, and getarg $$$$$
c
c   Asnag1 assigns logical unit lu to a direct access disk file
c   with mode "mode" and record length "len".  See dasign for 
c   details.  The n th argument is used as the model name.  If there 
c   is no n th argument and ib is non-blank, it is taken to be the 
c   model name.  If ib is blank, the user is prompted for the
c   model name using the character string in variable ia as the
c   prompt.  Programmed on 8 October 1980 by R. Buland.
c
      save
      logical log
      character*(*) ia,ib
      character cdum*30
c
c      if(iargc().lt.n) go to 1
c      call getarg(n,ib)
c      go to 2
c
c 1    if(ib.ne.' ') go to 2
c      call query(ia,log)
c      read(*,100)ib
c 100  format(a)
c
c 2    nb=index(ib,' ')-1
c      if(nb.le.0) nb=len(ib)
c      cdum=ib(1:nb)//'.hed'
      cdum='ak135.hed'
      call assign(lu,mode,cdum)
      return
      end
      subroutine assign(lu,mode,ia)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine assign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by assign are
c   assumed to be sequential.  Programmed on 3 December 1979 by
c   R. Buland.
c
      save
      character*(*) ia
      logical exst
c
      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.3) ns=3
      go to (1,2),nf
 1    go to (11,12,13),ns
 11   open(lu,file=ia,status='old',form='formatted')
      rewind lu
      return
 12   inquire(file=ia,exist=exst)
      if(exst) go to 11
 13   open(lu,file=ia,status='new',form='formatted')
      return
 2    go to (21,22,23),ns
 21   open(lu,file=ia,status='old',form='unformatted')
      rewind lu
      return
 22   inquire(file=ia,exist=exst)
      if(exst) go to 21
 23   open(lu,file=ia,status='new',form='unformatted')
      return
      end
      subroutine retrns(lu)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine retrns closes (disconnects) logical unit lu from the
c   calling program.  Programmed on 3 December 1979 by R. Buland.
c
      save
      close(unit=lu)
      return
      end
      subroutine query(ia,log)
c
c $$$$$ calls tnoua $$$$$
c
c   Subroutine query scans character string ia (up to 78 characters) for
c   a question mark or a colon.  It prints the string up to and
c   including the flag character plus two blanks with no newline on the
c   standard output.  If the flag was a question mark, query reads the
c   users response.  If the response is 'y' or 'yes', log is set to
c   true.  If the response is 'n' or 'no', log is set to false.  Any
c   other response causes the question to be repeated.  If the flag was
c   a colon, query simply returns allowing user input on the same line.
c   If there is no question mark or colon, the last non-blank character
c   is treated as if it were a colon.  If the string is null or all
c   blank, query prints an astrisk and returns.  Programmed on 3
c   December 1979 by R. Buland.
c
      save
      logical log
      character*(*) ia
      character*81 ib
      character*4 ans
      nn=len(ia)
      log=.true.
      ifl=1
      k=0
c   Scan ia for flag characters or end-of-string.
      do 1 i=1,nn
      ib(i:i)=ia(i:i)
      if(ib(i:i).eq.':') go to 7
      if(ib(i:i).eq.'?') go to 3
      if(ib(i:i).eq.'\0') go to 5
      if(ib(i:i).ne.' ') k=i
 1    continue
c   If we fell off the end of the string, branch if there were any non-
c   blank characters.
 5    if(k.gt.0) go to 6
c   Handle a null or all blank string.
      i=1
      ib(i:i)='*'
      go to 4
c   Handle a string with no question mark or colon but at least one
c   non-blank character.
 6    i=k
c   Append two blanks and print the string.
 7    i=i+2
      ib(i-1:i-1)=' '
      ib(i:i)=' '
c   Tnoua prints the first i characters of ib without a newline.
 4    call tnoua(ib,i)
      if(ifl.gt.0) return
c   If the string was a yes-no question read the response.
      read 102,ans
 102  format(a4)
      call uctolc(ans,-1)
c   If the response is yes log is already set properly.
      if(ans.eq.'y   '.or.ans.eq.'yes ') return
c   If the response is no set log to false.  Otherwise repeat the
c   question.
      if(ans.ne.'n   '.and.ans.ne.'no  ') go to 4
      log=.false.
      return
 3    ifl=-ifl
      go to 7
      end
      subroutine uctolc(ia,ifl)
c
c $$$$$ calls only library routines $$$$$
c
c   Subroutine uctolc converts alphabetic characters in string ia from
c   upper case to lower case.  If ifl < 0 all characters are converted.
c   Otherwise characters enclosed by single quotes are left unchanged.
c   Programmed on 21 January by R. Buland.  Calling sequence changed
c   on 11 December 1985 by R. Buland.
c
      character*(*) ia
      data nfl/1/
      if(ifl.lt.0) nfl=1
c   Scan the string.
      n=len(ia)
      do 1 i=1,n
      if(ifl.lt.0) go to 2
c   Look for single quotes.
      if(ia(i:i).eq.'''') nfl=-nfl
c   If we are in a quoted string skip the conversion.
      if(nfl.lt.0) go to 1
c   Do the conversion.
 2    ib=ichar(ia(i:i))
      if(ib.lt.65.or.ib.gt.90) go to 1
      ia(i:i)=char(ib+32)
 1    continue
      return
      end
      subroutine r4sort(n,rkey,iptr)
c
c $$$$$ calls no other routine $$$$$
c
c   R4sort sorts the n elements of array rkey so that rkey(i), 
c   i = 1, 2, 3, ..., n are in asending order.  R4sort is a trivial
c   modification of ACM algorithm 347:  "An efficient algorithm for
c   sorting with minimal storage" by R. C. Singleton.  Array rkey is
c   sorted in place in order n*alog2(n) operations.  Coded on
c   8 March 1979 by R. Buland.  Modified to handle real*4 data on
c   27 September 1983 by R. Buland.
c
      save
      dimension rkey(n),iptr(n),il(10),iu(10)
c   Note:  il and iu implement a stack containing the upper and
c   lower limits of subsequences to be sorted independently.  A
c   depth of k allows for n<=2**(k+1)-1.
      if(n.le.0) return
      do 1 i=1,n
 1    iptr(i)=i
      if(n.le.1) return
      r=.375
      m=1
      i=1
      j=n
c
c   The first section interchanges low element i, middle element ij,
c   and high element j so they are in order.
c
 5    if(i.ge.j) go to 70
 10   k=i
c   Use a floating point modification, r, of Singleton's bisection
c   strategy (suggested by R. Peto in his verification of the
c   algorithm for the ACM).
      if(r.gt..58984375) go to 11
      r=r+.0390625
      go to 12
 11   r=r-.21875
 12   ij=i+(j-i)*r
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 20
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 20   l=j
      if(rkey(iptr(j)).ge.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(j)
      iptr(j)=it
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 39   tmpkey=rkey(iptr(ij))
      go to 40
c
c   The second section continues this process.  K counts up from i and
c   l down from j.  Each time the k element is bigger than the ij
c   and the l element is less than the ij, then interchange the
c   k and l elements.  This continues until k and l meet.
c
 30   it=iptr(l)
      iptr(l)=iptr(k)
      iptr(k)=it
 40   l=l-1
      if(rkey(iptr(l)).gt.tmpkey) go to 40
 50   k=k+1
      if(rkey(iptr(k)).lt.tmpkey) go to 50
      if(k.le.l) go to 30
c
c   The third section considers the intervals i to l and k to j.  The
c   larger interval is saved on the stack (il and iu) and the smaller
c   is remapped into i and j for another shot at section one.
c
      if(l-i.le.j-k) go to 60
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 80
 60   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 80
c
c   The fourth section pops elements off the stack (into i and j).  If
c   necessary control is transfered back to section one for more
c   interchange sorting.  If not we fall through to section five.  Note
c   that the algorighm exits when the stack is empty.
c
 70   m=m-1
      if(m.eq.0) return
      i=il(m)
      j=iu(m)
 80   if(j-i.ge.11) go to 10
      if(i.eq.1) go to 5
      i=i-1
c
c   The fifth section is the end game.  Final sorting is accomplished
c   (within each subsequence popped off the stack) by rippling out
c   of order elements down to their proper positions.
c
 90   i=i+1
      if(i.eq.j) go to 70
      if(rkey(iptr(i)).le.rkey(iptr(i+1))) go to 90
      k=i
      kk=k+1
      ib=iptr(kk)
 100  iptr(kk)=iptr(k)
      kk=k
      k=k-1
      if(rkey(ib).lt.rkey(iptr(k))) go to 100
      iptr(kk)=ib
      go to 90
      end
      function iupcor(phnm,dtdd,xcor,tcor)
      save
      include 'ttlim.inc'
      character*(*) phnm
      character*8 phcd
      double precision us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,xu,
     1 px,xt,taut,coef,tauc,xc,tcoef,tp
      double precision x,dp,dps,ps,cn
      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      data oldep,jp,js/-1.,2*0/,cn/57.295779d0/
c
      iupcor=1
c     print *,'oldep odep',oldep,odep
      if(oldep.eq.odep) go to 1
      oldep=odep
c   Find the upgoing P branch.
c     print *,'mbr1 mbr2',mbr1,mbr2
      do 2 jp=mbr1,mbr2
c     print *,'jp phcd xbrn',jp,'  ',phcd(jp),xbrn(jp,1)
      if((phcd(jp).eq.'Pg'.or.phcd(jp).eq.'Pb'.or.phcd(jp).eq.'Pn'.or.
     1 phcd(jp).eq.'P').and.xbrn(jp,1).le.0d0) go to 3
 2    continue
      jp=0
c   Find the upgoing S branch.
 3    do 4 js=mbr1,mbr2
c     print *,'js phcd xbrn',js,'  ',phcd(js),xbrn(js,1)
      if((phcd(js).eq.'Sg'.or.phcd(js).eq.'Sb'.or.phcd(js).eq.'Sn'.or.
     1 phcd(js).eq.'S').and.xbrn(js,1).le.0d0) go to 1
 4    continue
      js=0
c
c1    print *,'jp js',jp,js
 1    if(phnm.ne.'P'.and.phnm.ne.'p') go to 5
      jb=jp
      if(jb)14,14,6
c
 5    if(phnm.ne.'S'.and.phnm.ne.'s') go to 13
      jb=js
      if(jb)14,14,6
c
 6    is=jndx(jb,1)+1
      ie=jndx(jb,2)
      ps=abs(dtdd)/dn
c     print *,'jb is ie dtdd dn ps',jb,is,ie,dtdd,dn,ps
      if(ps.lt.pt(is-1).or.ps.gt.pt(ie)) go to 13
      do 7 i=is,ie
c     print *,'i pt',i,pt(i)
      if(ps.le.pt(i)) go to 8
 7    continue
      go to 13
c
 8    j=i-1
      dp=pt(ie)-ps
      dps=dsqrt(dabs(dp))
      x=tau(2,j)+2d0*dp*tau(3,j)+1.5d0*dps*tau(4,j)
c     print *,'j pt dp dps x',j,pt(ie),dp,dps,x
      tcor=tn*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+dps*tau(4,j))+ps*x)
      xcor=cn*x
c     print *,'iupcor xcor tcor',iupcor,xcor,tcor
      return
c
 13   iupcor=-1
 14   xcor=0.
      tcor=0.
c     print *,'iupcor xcor tcor',iupcor,xcor,tcor
      return
      end
      subroutine brnset(nn,pcntl,prflg)
c
c   Brnset takes character array pcntl(nn) as a list of nn tokens to be
c   used to select desired generic branches.  Prflg(3) is the old
c   prnt(2) debug print flags in the first two elements plus a new print
c   flag which controls a branch selected summary from brnset.  Note that
c   the original two flags controlled a list of all tau interpolations
c   and a branch range summary respectively.  The original summary output
c   still goes to logical unit 10 (ttim1.lis) while the new output goes
c   to the standard output (so the caller can see what happened).  Each
c   token of pcntl may be either a generic branch name (e.g., P, PcP,
c   PKP, etc.) or a keyword (defined in the data statement for cmdcd
c   below) which translates to more than one generic branch names.  Note
c   that generic branch names and keywords may be mixed.  The keywords
c   'all' (for all branches) and 'query' (for an interactive token input
c   query mode) are also available.
c
      save
      parameter(ncmd=4,lcmd=16)
      include 'ttlim.inc'
      logical prflg(3),segmsk,prnt,fnd,all
      character*(*) pcntl(nn)
      character*8 phcd,segcd(jbrn),cmdcd(ncmd),cmdlst(lcmd),phtmp,
     1 phlst(jseg)
      double precision zs,pk,pu,pux,tauu,xu,px,xt,taut,coef,tauc,xc,
     1 tcoef,tp
      dimension nsgpt(jbrn),ncmpt(2,ncmd)
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
c   Segmsk is a logical array that actually implements the branch
c   editing in depset and depcor.
      common/prtflc/segmsk(jseg),prnt(2)
c
c   The keywords do the following:
c      P      gives P-up, P, Pdiff, PKP, and PKiKP
c      P+     gives P-up, P, Pdiff, PKP, PKiKP, PcP, pP, pPdiff, pPKP,
c             pPKiKP, sP, sPdiff, sPKP, and sPKiKP
c      S+     gives S-up, S, Sdiff, SKS, sS, sSdiff, sSKS, pS, pSdiff,
c             and pSKS
c      basic  gives P+ and S+ as well as ScP, SKP, PKKP, SKKP, PP, and
c             P'P'
c   Note that generic S gives S-up, Sdiff, and SKS already and so
c   doesn't require a keyword.
c
      data cmdcd/'P','P+','basic','S+'/
      data cmdlst/'P','PKiKP','PcP','pP','pPKiKP','sP','sPKiKP','ScP',
     1 'SKP','PKKP','SKKP','PP','S','ScS','sS','pS'/
      data ncmpt/1,2,1,7,1,13,13,16/
c
c   Take care of the print flags.
      prnt(1)=prflg(1)
      prnt(2)=prflg(2)
      if(prnt(1)) prnt(2)=.true.
c   Copy the token list into local storage.
      no=min0(nn,jseg)
      do 23 i=1,no
 23   phlst(i)=pcntl(i)
c   See if we are in query mode.
      if(no.gt.1.or.(phlst(1).ne.'query'.and.phlst(1).ne.'QUERY'))
     1 go to 1
c
c   In query mode, get the tokens interactively into local storage.
c
 22   print *,'Enter desired branch control list at the prompts:'
      no=0
 21   call query(' ',fnd)
      if(no.ge.jseg) go to 1
      no=no+1
      read 100,phlst(no)
 100  format(a)
c   Terminate the list of tokens with a blank entry.
      if(phlst(no).ne.' ') go to 21
      no=no-1
      if(no.gt.0) go to 1
c   If the first token is blank, help the user out.
      print *,'You must enter some branch control information!'
      print *,'     possibilities are:'
      print *,'          all'
      print 101,cmdcd
 101  format(11x,a)
      print *,'          or any generic phase name'
      go to 22
c
c   An 'all' keyword is easy as this is already the default.
 1    all=.false.
      if(no.eq.1.and.(phlst(1).eq.'all'.or.phlst(1).eq.'ALL'))
     1 all=.true.
      if(all.and..not.prflg(3)) return
c
c   Make one or two generic branch names for each segment.  For example,
c   the P segment will have the names P and PKP, the PcP segment will
c   have the name PcP, etc.
c
      kseg=0
      j=0
c   Loop over the segments.
      do 2 i=1,nseg
      if(.not.all) segmsk(i)=.false.
c   For each segment, loop over associated branches.
 9    j=j+1
      phtmp=phcd(j)
c   Turn the specific branch name into a generic name by stripping out
c   the crustal branch and core phase branch identifiers.
      do 3 l=2,8
 6    if(phtmp(l:l).eq.' ') go to 4
      if(phtmp(l:l).ne.'g'.and.phtmp(l:l).ne.'b'.and.phtmp(l:l).ne.'n')
     1 go to 5
      if(l.lt.8) phtmp(l:)=phtmp(l+1:)
      if(l.ge.8) phtmp(l:)=' '
      go to 6
 5    if(l.ge.8) go to 3
      if(phtmp(l:l+1).ne.'ab'.and.phtmp(l:l+1).ne.'ac'.and.
     1 phtmp(l:l+1).ne.'df') go to 3
      phtmp(l:)=' '
      go to 4
 3    continue
c4    print *,'j phcd phtmp =',j,' ',phcd(j),' ',phtmp
c
c   Make sure generic names are unique within a segment.
 4    if(kseg.lt.1) go to 7
      if(phtmp.eq.segcd(kseg)) go to 8
 7    kseg=kseg+1
      segcd(kseg)=phtmp
      nsgpt(kseg)=i
c     if(prflg(3)) print *,'kseg nsgpt segcd =',kseg,nsgpt(kseg),' ',
c    1 segcd(kseg)
 8    if(jidx(j).lt.indx(i,2)) go to 9
 2    continue
      if(all) go to 24
c
c   Interpret the tokens in terms of the generic branch names.
c
      do 10 i=1,no
c   Try for a keyword first.
      do 11 j=1,ncmd
      if(phlst(i).eq.cmdcd(j)) go to 12
 11   continue
c
c   If the token isn't a keyword, see if it is a generic branch name.
      fnd=.false.
      do 14 k=1,kseg
      if(phlst(i).ne.segcd(k)) go to 14
      fnd=.true.
      l=nsgpt(k)
      segmsk(l)=.true.
c     print *,'Brnset:  phase found - i k l segcd =',i,k,l,' ',
c    1 segcd(k)
 14   continue
c   If no matching entry is found, warn the caller.
      if(.not.fnd) print *,'Brnset:  phase ',phlst(i),' not found.'
      go to 10
c
c   If the token is a keyword, find the matching generic branch names.
 12   j1=ncmpt(1,j)
      j2=ncmpt(2,j)
      do 15 j=j1,j2
      do 15 k=1,kseg
      if(cmdlst(j).ne.segcd(k)) go to 15
      l=nsgpt(k)
      segmsk(l)=.true.
c     print *,'Brnset:  cmdlst found - j k l segcd =',j,k,l,' ',
c    1 segcd(k)
 15   continue
 10   continue
c
c   Make the caller a list of the generic branch names selected.
c
 24   if(.not.prflg(3)) return
      fnd=.false.
      j2=0
c   Loop over segments.
      do 16 i=1,nseg
      if(.not.segmsk(i)) go to 16
c   If selected, find the associated generic branch names.
      j2=j2+1
      do 17 j1=j2,kseg
      if(nsgpt(j1).eq.i) go to 18
 17   continue
      print *,'Brnset:  Segment pointer (',i,') missing?'
      go to 16
 18   do 19 j2=j1,kseg
      if(nsgpt(j2).ne.i) go to 20
 19   continue
      j2=kseg+1
c   Print the result.
 20   j2=j2-1
c     if(.not.fnd) print *,'Brnset:  the following phases have '//
c    1 'been selected -'
      fnd=.true.
c     print 102,i,(segcd(j),j=j1,j2)
c102  format(10x,i5,5(2x,a))
 16   continue
      return
      end
c
      block data

      include "ttlim.inc"

      logical segmsk,prnt
      common/prtflc/ segmsk(jseg),prnt(2)
      data segmsk,prnt/jseg*.true.,2*.false./ 

      double precision zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
c
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      data tauc,xc/jtsm*0d0,jxsm*0d0/

      integer ka
      real deplim
      double precision ua,taua
      common/pdec/ua(5,2),taua(5,2),deplim,ka
      data deplim,ka/1.1,4/ 

      end
c------------------------------------------------------------------------
c
c	nn2d_setup - Performs all setup procedures for natural neighbour 
c		     interpolation routine nn2D.
c
c	Input:
c	        np			number of nodes
c	        nt_max			maximum number of triangles
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		np_max			maximum number of nodes 
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~20 in calling program)
c		nmax			maximum sum of the number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c		points(2,np)		array of node co-ordinates
c	        dmode			Delaunay calculation mode (integer)
c	        nmode			NN setup mode
c	        clockwise		logical for the vertices input order  
c               data(np)                data values at each node
c		nnn			Integer work array used by build_nv
c		nnlist			Integer work array used by build_nv
c		ntwork		        Integer work array used by build_nv
c		nohalt_hull		determines error response in routine
c					calculate_hulltriangles
c		eps 			tolerance parameter used by delaun
c					(see delaun for details)
c		vis_tlist		Integer work array used by delaun
c		vis_elist		Integer work array used by delaun
c		add_tlist		Integer work array used by delaun
c		nv_max			size of delaun work arrays vis_tlist
c					vis_elist, & add_tlist (passed to 
c					delaun for error checking)
c
c	Output:
c	        nt			number of triangles
c	        vertices(3,nt)		array of triangle vertices 
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        nh			number of triangles with an edge
c					on the convex hull
c		hulltriangles(nh)	array of triangles with an edge 
c					on the convex hull
c		loc			an initial guess triangle for point 
c					location routine `triloc' used by nn2D
c					(set somewhere near the centre)
c
c	Operation modes:
c
c		 The setup routine will perform different tasks depending
c		 on the input parameters dmode and nmode (see table below).
c		 Depending on the modes used some work arrays may be set 
c		 to size 1 to save memory. The "Memory Savings" column in the
c		 table below shows the dimension statement that may
c		 be used in the calling program if the routine is ONLY EVER 
c		 CALLED IN THE CORRESPONDING MODE. 
c
c		 PARAMETERS	ACTION	 		MEMORY SAVINGS
c
c		 nmode = 1	Delaunay only 		real*8 centres(3,1)
c							integer hulltriangles(1)
c		 nmode = 0	Delaunay + nn setup	
c		 nmode = -1	nn setup only		
c
c		 dmode > 0	Delaunay read in from   integer vis_tlist(1)
c				logical unit dmode.	integer vis_elist(1) 
c							integer add_tlist(1) 
c		 dmode = 0      Qhull used		Same as dmode > 0.
c		 dmode = -1	Delaun + X-sort         integer nnn(1)
c							integer nnlist(1)
c							integer ntwork(1)
c		 dmode = -2	Delaun + no sort        Same as dmode=-1
c
c		 dmode = 0 & nmode=1			integer neighbour(3,1)
c		 
c		 A call with nmode = -1 can only be made after a call 
c		 with nmode = 1.
c
c	Comments:
c
c		 If the arrays are used then they should be dimensioned
c		 in the calling program in the following way:
c
c		 real*8		points(2,np_max)
c		 real*8		centres(3,nt_max)
c		 integer	vertices(3,nt_max)
c		 integer	neighbour(3,nt_max)
c		 integer	hulltriangles(nh_max)
c		 
c		 Note: nh_max can usually be dimensioned much less than nt_max
c		 because hulltriangles, stores in a compact form, all 
c		 triangles with an edge on the convex hull. Except for
c		 very irregular point distributions nh << nt. If nh is 
c		 determined to be > nh_max then an error is reported and the
c		 program is halted (unless nohalt parameter is set to 1). 
c		 The array hulltriangles is only used by nn2Do see routine 
c		 calculate_hulltriangles. If nh_max = 1 then hulltriangles 
c		 is not calculated.
c
c		 The initial guess triangle 'loc' is set in nn_setup but
c		 at each call it will be set to the triangle found during
c		 the previous call to nn2D. The user may modify its value
c		 if the input point (x,y) is known to be in, or near, a
c		 particular triangle.
c
c		 If dmode > 0 the the deluanay tessellation is read in from
c		 logical unit `dmode' instead of being calculated internally 
c		 This can be useful if qhullf fails because
c		 of precision errors. The Deluanay may be determined
c		 externally to this program using a double precision version
c		 or another algorithm, e.g. Fortune's sweepline method.
c
c		 If 50 > dmode > 0 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ZERO, 
c		 which is the standard output format of codes qhull 
c		 (quickhull method) and voronoi (sweepline method).
c		 If clockwise = .true. (.false.) then the vertices are assumed 
c		 to be in clockwise (anti-clockwise) order. Note program
c		 qhull outputs vertices in anti-clockwise order while 
c		 voronoi in clockwise order. The internal format is 
c		 anti-clockwise and nodes numbered from ONE.
c
c		 If dmode => 50 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ONE,
c		 which is the output format of program del (using delaun). 
c
c		 Three other work arrays are produced as a `by product'
c		 of the routine build_nv which calculates the neighbour
c		 array. These must be dimensioned in the calling program in 
c		 the following way (unless delaun is used for calculating the
c		 Delaunay because it already determines the neighbour array)
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntwork(nmax) : dummy work array
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 Each of these are useful lists that describe features of
c		 the Voronoi diagram. Both nnlist and ntwork are stored in
c		 a compact format to avoid zeros. They are only used 
c		 in the setup routine and the memory may be freed once
c		 initialization is completed.
c		 
c
c		 Calls are made to: qhullf, ccentres, build_nv and 
c				    calculate_hulltriangles, delaun.
c		 
c					M. Sambridge, RSES, April 1994.
c  					        (Last modified 10/4/96)
c
c------------------------------------------------------------------------
c
	Subroutine nn2d_setup
     &             (np,nt_max,nh_max,np_max,nnpn_max,nmax,
     &              points,dmode,nmode,clockwise,data,nt,vertices,
     &              centres,neighbour,nh,hulltriangles,nohalt_hull,
     &              loc,nnn,nnlist,ntwork,
     &              eps,nv_max,vis_tlist,vis_elist,add_tlist,
     &		    lt_work,ln_work)

	real*8		points(2,*)
	real*8		centres(3,*)
        real*8          data(*)
	real*8		eps
	integer		vertices(3,*)
	integer		neighbour(3,*)
        integer		hulltriangles(*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntwork(*)
        integer         vis_tlist(*)
        integer         vis_elist(*)
        integer         add_tlist(*)
        integer         dmode,nmode
	logical*1	lt_work(*)
	logical*1	ln_work(*)
	logical		nnwrite
	logical		clockwise
	logical		ldummy
        logical         timing

        common/nnswitches/nnwrite,lud

        common/timing/timing,t_loc,t_int,t_setup

        if(timing)a = cputime(t1,t2)

        if(nmode.eq.1.or.nmode.eq.0)then

           if(dmode.eq.0)then
c                                       calculate Delaunay using qhull 
 
              call qhullf(np,2,2,nt_max,0,points,nt,vertices)

           else if(dmode.eq.-1.or.dmode.eq.-2)then

c					sort the points in ascending x order
c					and rearrange data points also
	      if(dmode.eq.-1)then

c		  write(*,*)' X sort in progress'
	          call hpsort_d(np,1,points,data)
c		  write(*,*)' X sort done'
 		  write(*,*)' Input points have been sorted in order '
 		  write(*,*)' of first co-ordinate'

	      end if

c                                       calculate Delaunay using delaun 
 
              call delaun (points,np,neighbour,vertices,nt,nt_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     0,ldummy,0,0,0)

	   else
c                                       read in Delaunay vertices

              nt = 0
              i1 = 1
              i2 = 2
              if(clockwise)then
                 i1 = 2
                 i2 = 1
              end if
              read(dmode,*)
  1           read(dmode,*,end=3,err=2)
     &        vertices(i1,nt+1),vertices(i2,nt+1),vertices(3,nt+1)
              nt = nt + 1
              if(nt.ge.nt_max)then
                 write(*,*) 'Error in nn_setup: too many triangles'
                 write(*,*) 'Remedy: increase size of parameter nt_max'
                 write(*,*) '        in calling program.'
                 stop 
              end if
              go to 1
  2           write(*,*)
     &        'Error in nn_setup: read error in Delaunay input file'
              stop
  3           continue
     
           end if
 
c					adjust array vertices to
c					range from nodes 1 to np
c
	   if(dmode.ge.0.and.dmode.lt.50)then
	      do 5 i = 1,nt
	         vertices(1,i) = vertices(1,i) + 1
	         vertices(2,i) = vertices(2,i) + 1
	         vertices(3,i) = vertices(3,i) + 1
 5            continue
	   end if

	end if
c
c					Perform set up for nn interpolation
c
        if(nmode.eq.0.or.nmode.eq.-1)then

c					set initial guess for 
c					triangle location procedure
           loc = nt/2

c                                       Calculate Circumcentres

           call ccentres(points,vertices,nt,centres)

c                                       Build neighbour matrix
c					(if not already built)

           if(dmode.ge.0)then
              call build_nv
     &        (np,vertices,nt,np_max,nmax,
     &         neighbour,nnn,nnlist,ntwork)
           end if

c					calculate hulltriangles

           if(nh_max.gt.1) call calculate_hulltriangles
     &     (neighbour,nt,nh_max,hulltriangles,nh,nohalt_hull)

c					initialize logical work arrays 
           do i=1,nt
              lt_work(i) = .false.
           end do
           do i=1,np
              ln_work(i) = .false.
           end do

	end if

        if(timing)then
           a = cputime(t1,t2)
           t_setup = t_setup + t1
        end if

	return
	end
c
c------------------------------------------------------------------------
c

c------------------------------------------------------------------------
c
c	Heapsort - Sorts an array into ascending order.
c		   Modified from numerical recipes 2, to use a 
c		   two dimensional double precision array.
c
c		   Sorting is done on index ka (ka=1 or 2).
c
c------------------------------------------------------------------------
c
c
      SUBROUTINE hpsort_d(n,ka,ra,da)
      INTEGER n
c     REAL ra(n)
      REAL*8 ra(2,n)
      REAL*8 da(n)
      INTEGER i,ir,j,l,ka
c     REAL rra
      REAL*8 rra(2)
      REAL*8 dda
      kb = 1
      if(ka.eq.1)kb = 2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra(ka)=ra(ka,l)
          rra(kb)=ra(kb,l)
	  dda = da(l)
        else
          rra(ka)=ra(ka,ir)
          rra(kb)=ra(kb,ir)
	  dda = da(ir)
          ra(ka,ir)=ra(ka,1)
          ra(kb,ir)=ra(kb,1)
          da(ir)=da(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(ka,1)=rra(ka)
            ra(kb,1)=rra(kb)
            da(1)=dda
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(ka,j).lt.ra(ka,j+1))j=j+1
          endif
          if(rra(ka).lt.ra(ka,j))then
            ra(ka,i)=ra(ka,j)
            ra(kb,i)=ra(kb,j)
            da(i)=da(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(ka,i)=rra(ka)
        ra(kb,i)=rra(kb)
        da(i)=dda
      goto 10
      END
c
c------------------------------------------------------------------------
c
c	pangle - pseudo angle routine 
c
c		 returns a number between 0 and 360 which is NOT the
c		 angle made by the line from p1 to p2 with the horizontal
c		 but which has the same order properties as that angle,
c		 i.e. has the same order of angles as arctan dy/dx.
c	 	 This function involves only simple products and quotients.
c
c		 From Sedgewick (1990) `Algorithms in C' (Addison Wesley)
c
c						M. Sambridge 1996.
c
c------------------------------------------------------------------------
c
	Function pangle(p1,p2)

	real*8		p1(2)
	real*8		p2(2)

	dx = p2(1) - p1(1)
        ax = abs(dx)
	dy = p2(2) - p1(2)
        ay = abs(dy)
        t = 0.
        a = ax+ay
	if(a.ne.0.)then
           t = dy/a
        end if
        if(dx.lt.0.)then
          t = 2-t
        else if(dy.lt.0.)then 
          t = 4+t
        end if
        pangle = t*90

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum - calculates circum-centre of three points
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c	Comments:
c
c		 Solves 3x3 linear system of equations.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum(pa,pb,pc,centre)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_d - calculates circum-centre of three points and 
c	           derivatives of circum-centre with respect to 
c		   co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_d(pa,pb,pc,centre,vx,vy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2)
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2

c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_dd - calculates circum-centre of three points and 
c		    1st and 2nd derivatives of circum-centre with 
c		    respect to co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_dd(pa,pb,pc,centre,vx,vy,vxx,vyy,vxy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2),vxx(2),vyy(2),vxy(2)
c
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2
c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2
c						Second derivatives
	   f11 = 2*vx(1) - 1.
	   f22 = 2*vy(2) - 1.
	   f12 = vy(1) + vx(2)

           vxx(1) = f11*denum1
           vxx(2) = f11*denum2
           vyy(1) = f22*denum1
           vyy(2) = f22*denum2
           vxy(1) = f12*denum1
           vxy(2) = f12*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c						heapsort modified to 
c						sort two arrays  
c
c------------------------------------------------------------------------
c
      SUBROUTINE hpsort_two(n,ra,rb)
      INTEGER n
      REAL ra(n)
      REAL*8 rb(2,n)
      INTEGER i,ir,j,l
      REAL rra
      REAL*8 rrb1,rrb2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          rrb1=rb(1,l)
          rrb2=rb(2,l)
        else
          rra=ra(ir)
          rrb1=rb(1,ir)
          rrb2=rb(2,ir)
          ra(ir)=ra(1)
          rb(1,ir)=rb(1,1)
          rb(2,ir)=rb(2,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            rb(1,1)=rrb1
            rb(2,1)=rrb2
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            rb(1,i)=rb(1,j)
            rb(2,i)=rb(2,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
        rb(1,i)=rrb1
        rb(2,i)=rrb2
      goto 10
      END
c
c------------------------------------------------------------------------
c
c					Numerical Recipes routine index
c
c------------------------------------------------------------------------
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
c
c------------------------------------------------------------------------
c
c	second_v_area - calculates the area of a second-order Voronoi 
c			cell using an un-ordered list of vertices and a
c		        a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		n			number of vertices
c
c	Output:
c		area			area of polygon X 2
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Comments:
c		 The vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area(x,n,p,a,area)

      real*8	x(2)
      real*8	p(2,*)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
      area = abs(area)
      
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_d - calculates the area of a second-order Voronoi 
c			  cell using an un-ordered list of vertices and a
c		          a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first two vertices 
c					of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx,
c					      df(2) = df/dy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_d(x,n,p,dp,a,area,df)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	df(2)
      real*4    theta
      real*4    a(*)
c
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
      end if

      area = abs(area)
       
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_dd - calculates the area of a second-order Voronoi 
c			   cell using an un-ordered list of vertices and a
c		           a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first 
c					two vertices of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx
c					      df(2) = df/dy
c		ddf(3)			second derivatives of area X 2
c					      ddf(1) = d2f/dxx
c					      ddf(2) = d2f/dyy
c					      ddf(3) = d2f/dxy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two and xplot routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_dd(x,n,p,dp,ddp,a,area,df,ddf)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	ddp(6,2)
      real*8	df(2),ddf(3)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
        a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

c						calculate 2nd derivatives
     
      ddf(1) =   ddp(1,1)*d222n + ddp(2,1)*d1n12 
     &         + ddp(1,2)*d2321 + ddp(2,2)*d1113
     &         + 2.*(dp(1,1)*dp(2,2) - dp(2,1)*dp(1,2))
      ddf(2) =   ddp(3,1)*d222n + ddp(4,1)*d1n12 
     &         + ddp(3,2)*d2321 + ddp(4,2)*d1113
     &         + 2.*(dp(3,1)*dp(4,2) - dp(4,1)*dp(3,2))
      ddf(3) =   ddp(5,1)*d222n + ddp(6,1)*d1n12 
     &         + ddp(5,2)*d2321 + ddp(6,2)*d1113
     &         + dp(1,1)*dp(4,2) - dp(2,1)*dp(3,2)
     &         + dp(3,1)*dp(2,2) - dp(4,1)*dp(1,2)

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
        ddf(1) = -ddf(1)
        ddf(2) = -ddf(2)
        ddf(3) = -ddf(3)
      end if

      area = abs(area)
       
      return
      end
c
c------------------------------------------------------------------------
c
c	delaun - calculates delaunay triangulation incrementally 
c	 	 for a set of points in 2-D using a variation of
c		 Lawson's algorithm.
c
c	Input:
c		points(2,np)		array of node co-ordinates
c		num			number of nodes to be used
c               vis_tlist(nv_max)       List of triangles visible from 
c					current point.
c               vis_elist(nv_max)       List of edges visible from 
c					current point.
c               add_tlist(nv_max)       work array used by routine addpoint
c               eps                     distance from an interface for a
c                                       a point to be considered on an
c                                       interface (real*8). Prevents zero
c                                       area triangles resulting from rounding
c                                       error when nodes are co-linear.
c		nv_max			size of work arrays
c		mode			(=0,1,2,3) operation mode (see below)
c		inactive(np)		logical array. If mode=1 then the i-th
c					node is ignored if active(i) = .true.
c		nfirst			If mode=3 then nfirst is the first
c					node added to an existing triangulation
c		itstart			If mode=3 then itstart is a first
c					guess triangle containing node first node
c		subset(np)		logical array. If mode=2 then only
c					the nodes (subset(i),i=1,num) are used.
c
c	Output:
c               v(3,*)           	array of triangle vertices
c               numtri                  number of triangles in current
c                                       triangulation.
c               e(3,*)                  adjacency matrix of neighbouring
c                                       triangles. e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2) 
c					in triangle j, stored counterclockwise 
c					about j.  
c                                       (This is the `opposite' definition)
c
c	Comments:
c
c       This routine calculates the Delaunay triangulation of a set of nodes 
c	using a variation of Lawson's method.  Each node is added sequentially 
c	and the Delaunay triangulation is updated. If the new node is inside 
c	the convex hull of the existing triangulation then the standard Lawson 
c	method is used. If it is outside then the list of triangle edges 
c	which are visible from the new point is calculated using routine 
c	visiblelist and each of these is used as the start of the swapping 
c	routine addpoint.
c
c	Four different operation modes are allowed.
c
c	MODE = 0:
c	The `standard' mode. All nodes from 1 to num are included. The arrays 
c	`subset' and `inactive' are unused and may be set to dummy variables
c       (saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 1:
c	All nodes from 1 to num are included except those for which
c	inactive(i) is set to true. The array `subset' is unused and may
c	be set to a dummy variable (saving memory). The variables nfirst 
c	and itstart are also unused.
c
c	MODE = 2:
c	Only nodes from subset(1) to subset(num) are included. The array
c	`inactive' is unused and may be set to a dummy variable
c	(saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 3:
c	Used to add nodes from nfirst to num to an existing triangulation.
c	Nodes for which inactive(i) is set to true are ignored.
c	The array `subset' is unused and may be set to a dummy variable.
c
c       The performance may be sensitive to the order in which the nodes are
c       added so these can be sorted before calling this routine if desired.
c
c	This routine was converted to use the `opposite' definition of
c	the adjacency matrix on 30/1/96.
c
c
c	Calls are made to triloc_del,visiblelist,insert_point,addpoint.
c
c					         M. Sambridge, Dec. 1994.
c					Modified by J. Braun, Sept. 1995.
c					(last change 30/1/96: multiple modes,
c					 and uses opposite definition of
c					 adjacency matrix)
c
c------------------------------------------------------------------------
c
	Subroutine delaun (points,num,e,v,numtri,numtri_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     mode,inactive,nfirst,itstart,subset)

	real*8		points(2,*)
	real*8		x,y
        real*8          eps,del1,del2,del
 	integer		vis_tlist(*),vis_elist(*),add_tlist(*)
 	integer		v(3,*)
 	integer		e(3,*)
 	integer		subset(*)
	integer		t,p,ccw
        logical         out
	logical		newpoint
	logical		inactive(*)

        if (mode.eq.0.or.mode.eq.1.or.mode.eq.2) then

c					We are calculating Delaunay 
c					of all input points or a
c					subset of all input points

c					find first two active nodes
           if(mode.eq.0)then
              i1=1 
              i2=2 
              nodestart=3
           else if(mode.eq.1)then
              i1=0
              i2=0
              do i=1,num
                 if (i2.ne.0) goto 2222
                 if (i1.ne.0.and..not.inactive(i)) i2=i
                 if (i1.eq.0.and..not.inactive(i)) i1=i
              end do
 2222         continue
              nodestart = i2+1
           else if(mode.eq.2)then
              i1 = subset(1)
              i2 = subset(2)
              nodestart = 3
           end if
c                                       Find three non-colinear points
c                                       to form the first triangle 
           v(1,1) = i1
           v(2,1) = i2
           do 10 j=nodestart,num
              i = j
              if(mode.eq.2)then
                  i = subset(j)
              else if(mode.eq.1.and.inactive(i))then
                  go to 10
              end if
              istart=i
	      del1 = (points(2,i1)-points(2,i))
     &              *(points(1,i2)-points(1,i))
	      del2 = (points(1,i1)-points(1,i))
     &              *(points(2,i2)-points(2,i))
              del = del1-del2
              if(dabs(del).gt.eps) goto 11111
 10        continue
           stop 'all input data are in a line...'
11111      v(3,1) = istart

c					Initialize adjacency matrix
 	   e(1,1) = 0
 	   e(2,1) = 0
 	   e(3,1) = 0
c					Ensure initial triangle 
c					is in ccw order
c					
 	   if(ccw(points(1,v(1,1)),
     &            points(1,v(2,1)),
     &            points(1,v(3,1)),k).eq.-1)then
                  itemp = v(1,1)
                  v(1,1) = v(2,1)
                  v(2,1) = itemp
c	          write(*,*)' initial triangle was cw'
 	    end if
	
c					Initialize variables
 	    numtri = 1
	    t = 1

        else if (mode.eq.3) then
c					We are adding nodes to an 
c					existing triangulation
c					Perform initialization
           nodestart=nfirst
           t = itstart
           istart = 0
           if(t.le.0.or.t.gt.numtri)t=1

        end if 
c					Incrementally update the 
c					Delaunay triangulation

 	do 100 j=nodestart,num


           p = j
           if(mode.eq.1.or.mode.eq.3)then
             if (inactive(j)) goto 100
           else if(mode.eq.2)then
	     p = subset(j)
           end if
             
           if(p.eq.istart)go to 100

	   x = points(1,p)
	   y = points(2,p)

c					locate triangle 
c					containing current node

	   call triloc_del(x,y,points,v,e,t,eps,out,ipos,iface)


	   if(out)then

c					point is outside of convex hull, 
c					so find list of edges that are 
c					visible from current point 

	      call visiblelist(points,e,v,x,y,t,ipos,eps,
     &                         vis_tlist,vis_elist,nvis)

c					for each visible edge 
c					start swapping algorithm

              newpoint = .true.

	      if(nvis.gt.nv_max)then
                 write(*,*)' Error in subroutine delaun:'
                 write(*,*)' Too many visible triangles
     &                       from current point'
                 write(*,*)' Remedy: increase size of parameter nv_max'
                 write(*,*)'         in calling program'
                 write(*,*)'         Number of visible triangles '
                 write(*,*)'         for this point             =',nvis
                 write(*,*)'       Current value of nv_max    =',nv_max
                 stop
              end if

	      do 60 i=1,nvis
                 t = vis_tlist(i)
                 ipos = vis_elist(i)
                 jpos = mod(vis_elist(i),3)+1
c	         write(6,*)' visible t =',t,' node',v(ipos,t),v(jpos,t)
 	         call addpoint
     &           (points,e,v,p,t,ipos,numtri,newpoint,add_tlist) 
                 newpoint = .false.
 60           continue

	   else
 
c	      write(6,*)' point located in triangle',t

c					add node to inside of convex hull
c					using swapping algorithm

  	      call insertpoint(points,e,v,p,t,numtri,iface) 

           end if

           if (numtri.gt.numtri_max) then
              write (*,*) 'Error in subroutine delaun:'
              write (*,*) 'Too many triangles'
              write(*,*)' Remedy: increase size of parameter numtri_max'
              write(*,*)'         in calling program'
              write(*,*)'         Number of triangles '
              write(*,*)'         for this point             =',numtri
              write(*,*)'         Current value of numtri_max    =',
     &                  numtrimax_max
              stop
           endif

 100    continue

	return
	end
c------------------------------------------------------------------------
c
c	Subroutine visiblelist - calculates all sides of triangles 
c		                 visible from the point (x,y), which is
c			         outside of the convex hull.
c			
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		x,y			Co-ordinates of test point p
c		t			Index of any triangle on hull 
c					that is visible from point p.
c					(Usually given by routine Triloc_del.)
c		tpos			Position of edge in triangle t
c					(using Sloan's adjacency convention)
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c		nvis			Number of triangles visible from point p
c		vis_tlist		List of triangles visible from p
c		vis_elist		List of edges visible from p
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c		 This routine was converted from using Sloan's definition of
c		 the adjacency matrix to the `opposite' definition on 30/1/96.
c
c	Calls routine visible.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c		
	Subroutine visiblelist
     &             (points,neighbour,vertices,x,y,t,tpos,eps,
     &              vis_tlist,vis_elist,nvis)

	real*8		points(2,*)
	real*8		x,y
	real*8		eps
 	integer		vertices(3,*)
 	integer		neighbour(3,*)
 	integer		vis_tlist(*),vis_elist(*)
 	integer		t,tpos,pos,t1,t2,edg,tnew
	logical		visible
	logical		special
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	nvis = 1
	vis_tlist(1) = t
	vis_elist(1) = tpos
    	inode = vertices(tpos,t)
cd      write(6,100)t,inode,vertices(mod(tpos,3)+1,t)
    	pos = c1(tpos)
      	jnode = vertices(pos,t)
    	t1 = neighbour(c2(pos),t)
        special = .false.
	if(t1.eq.0)then
	  t1 = t
          tnew = 0
          special = .true.
	end if

  5     continue
        if(.not.special)then
           pos = edg(t1,jnode,vertices)
           tnew = neighbour(c2(pos),t1)
cd         write(6,*)' tnew =',tnew,' t1',t1,' jnode',jnode
        end if
        special = .false.
        if(tnew.eq.0)then
  6        continue
  	   if(visible(x,y,points,vertices,t1,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t1
	      vis_elist(nvis) = pos
cd            write(6,100)t1,jnode,vertices(mod(pos,3)+1,t1)
	   else
cd	      write(6,200)t1,jnode,vertices(mod(pos,3)+1,t1)
              go to 10
	   end if
           pos = c1(pos)
	   jnode = vertices(pos,t1)
    	   tnew = neighbour(c2(pos),t1)
	   if(tnew.eq.0) go to 6
           t1 = tnew
	   go to 5 
        else
cd	   write(6,300)t1,jnode,vertices(mod(pos,3)+1,t1)
	   t1 = tnew
	   go to 5 
	end if

  10	jnode = inode
    	pos = c2(tpos)
    	t2 = neighbour(c2(pos),t)
        special = .false.
	if(t2.eq.0)then
	  t2 = t
          tnew = 0
          special = .true.
	end if

  15    continue
        if(.not.special)then
           pos = c2(edg(t2,jnode,vertices))
           tnew = neighbour(c2(pos),t2)
        end if
        special = .false.
        if(tnew.eq.0)then
  16       continue
  	   if(visible(x,y,points,vertices,t2,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t2
	      vis_elist(nvis) = pos
cd            write(6,100)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
	   else
cd	      write(6,200)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
              go to 20
	   end if
	   jnode = vertices(pos,t2)
    	   pos = c2(pos)
    	   tnew = neighbour(c2(pos),t2)
	   if(tnew.eq.0)go to 16
	   t2 = tnew
	   go to 15
        else
cd	   write(6,300)t2,jnode,vertices(mod(pos,3)+1,t2)
	   t2 = tnew
	   go to 15
	end if

 20	continue
	      
c 100    format
c     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is visible')
c 200    format
c     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not visible')
c 300    format
c     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not on convex hull')

	return
	end
c
c------------------------------------------------------------------------
c
c	Function visible - determines whether the triangle t is visible
c		           from the point p on edge tpos.
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        t			Triangle to be tested
c	        tpos			Edge to be tested in triangle t 
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        visible			Logical: = true if edge is visible  
c	                                         = false if edge is not visible
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c	Calls no other routines.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Function visible(x,y,points,vertices,t,tpos,eps)

	real*8		points(2,*)
	real*8		del1,del2
	real*8		x,y
 	integer		vertices(3,*)
 	integer		t,tpos
	logical		visible
	real*8		eps,del
	integer		c1(3)
        save            c1
	data		c1/2,3,1/

        j = c1(tpos)
c						test edge tpos in triangle t
        i1 = vertices(tpos,t)
        i2 = vertices(j,t)
        del1 = (points(2,i1)-y)*(points(1,i2)-x)
        del2 = (points(1,i1)-x)*(points(2,i2)-y)
	del = del1-del2
        if(del.gt.eps)then
           visible = .true.
	else
           visible = .false.
	end if

	return
	end
c
c------------------------------------------------------------------------
c
c	addpoint - inserts a point into an existing delaunay triangulation 
c		   when point is outside of triangulation (but attached to
c		   triangle t) using the stacking procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle on convex hull visible 
c					from input point 
c		numtri			number of triangles in current
c					triangulation.
c		tpos			position of start node in triangle t
c		tri			list of triangles visible from point p
c		newpoint		logical = true if t is the first
c					triangle on the hull visible from p
c
c	Output:
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	When newpoint = false then there are multiple triangles
c	from the new point to the convex hull, and addpoint must
c	be called once for each attached triangle. In this case
c	the initialization of the adjacency list includes the
c	neighouring triangles already processed by addpoint, i.e.
c	those from point p to the hull.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine addpoint (points,e,v,p,t,tpos,numtri,newpoint,tri) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t,tpos,ccw
	logical		swap
	integer		tri(*)
	logical		newpoint
        save 		ip,lp
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	if(newpoint)then
	   ip = 0
	   lp = 0
        end if

c			Add new node to existing triangulation

c					create new triangle
        numtri = numtri + 1
        v1 = v(tpos,t)
        v2 = v(c1(tpos),t)
        if(ccw(points(1,v1),
     &         points(1,v2),
     &         points(1,p),k).eq.-1)then
               itemp = v1
               v1 = v2
               v2 = itemp
        end if
        v(1,numtri) = p
        v(2,numtri) = v1
        v(3,numtri) = v2

c					initialize adjacency list including
c					neighbouring triangles attached
c					from the point to the hull.

        e(c2(1),numtri) = 0
        e(c2(2),numtri) = t
        e(c2(3),numtri) = 0

c				
        if(.not.newpoint)then
           do 10 j=1,lp
              k = tri(j)
              if(v(2,k).eq.v1)then
c                write(6,*)' v1 match with node 2'
c                write(6,*)' current triangle',numtri,' new',k
c                write(6,*)' nodes:',v(1,k),v(2,k),v(3,k) 
c                write(6,*)' e mat:',e(c2(1),k),e(c2(2),k),e(c2(3),k) 
                 e(c2(1),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v1)then
                 e(c2(1),numtri) = k
                 e(c2(3),k) = numtri
              end if
              if(v(2,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(3),k) = numtri
              end if
 10        continue
        end if

c
c					initialize stack

 	call stackinit

c                                       update adjacency list
c                                       for triangle on old boundary
        e(c2(tpos),t) = numtri


c					add new triangle on stack

	call push(numtri)

c					loop while stack is not empty

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl) 
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = l
	   end if

	   if(b.ne.0)then
              call push(r)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = r
           end if
	   if(c.ne.0) e(edg(c,l,e),c)=r

	else

c					record triangle attached to p
	   ip = ip + 1
           tri(ip) = l

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
        call stackflush()

	lp = ip

c       write(6,*)' Number of triangles attached to last point',ip
c1	write(6,*)(tri(i),i=1,ip)
c	write(6,*)' triangles attached to last point on hull'
c       do 100 i=1,ip
c          it=tri(i)
c          do 101 k=1,3
c             l=mod(k,3)+1
c             if(e(k,it).eq.0)then
c                write(6,*)' t',it,' edge ',v(k,it),v(l,it)
c             end if
c 101      continue 
c 100   continue 
c

	return
	end
c
c------------------------------------------------------------------------
c
c	insertpoint - inserts a point into an existing delaunay triangulation 
c		      (when new point is inside triangle t) using the stacking 
c		      procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle containing input point
c		numtri			number of triangles in current
c					triangulation.
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Output:
c
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The new point is assumed to be inside the convex hull of the
c	existing triangulation.
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine insertpoint (points,e,v,p,t,numtri,iface) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t
	logical		swap
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/


c					add new node to existing triangulation
        if(iface.eq.0)then
	   a = e(c2(1),t)
	   b = e(c2(2),t)
	   c = e(c2(3),t)
	   v1 = v(1,t)
	   v2 = v(2,t)
	   v3 = v(3,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+2	
	   e(c2(2),t) = a	
	   e(c2(3),t) = numtri+1	
        

c					create new triangles
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v2
	   v(3,numtri)=v3
	   e(c2(1),numtri)=t
	   e(c2(2),numtri)=b
	   e(c2(3),numtri)=numtri+1
	   numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=numtri-1
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        else
           j = iface
           k = c1(j)
           i = c1(k)
	   a = e(c2(i),t)
	   b = e(c2(j),t)
	   c = e(c2(k),t)
	   v1 = v(i,t)
	   v2 = v(j,t)
	   v3 = v(k,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+1
	   e(c2(2),t) = a	
	   e(c2(3),t) = 0 	

c					create new triangle
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=0
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        end if
  
c
c					initialize stack

 	call stackinit

c					add new triangles on stack
c					and update adjacency list

	if(a.ne.0)call push(t)

	if(b.ne.0)then
          e(edg(b,t,e),b)=numtri-1
          call push(numtri-1)
        end if

	if(c.ne.0)then
          e(edg(c,t,e),c)=numtri
          call push(numtri)
        end if
c					loop while stack is not empty

        if(a.eq.0.and.b.eq.0.and.c.eq.0)go to 100

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl)
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   end if
	   if(b.ne.0) call push(r)
	   if(c.ne.0) e(edg(c,l,e),c)=r

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
 100    continue
        call stackflush()

	return
	end
c
c------------------------------------------------------------------------
c
c	Function edg - finds edge in triangle l which is adjacent 
c		       to triangle k.
c
c		       (From Sloan 1987)
c
c------------------------------------------------------------------------
c
	Function edg(l,k,e)
c
	integer		l,k,i,e(3,*),edg
c
	do 10 i=1,3
	   if(e(i,l).eq.k)then
              edg = i
              return
           end if
 10     continue

	write(*,*)' ***Error in function edg***'
	write(*,*)' ***Triangles not adjacent***'
	write(*,*)' triangle = ',l,' looking for triangle',k

	stop
	end
c
c------------------------------------------------------------------------
c
c	logical function swap - checks to see if point p lies 
c			        inside circumcircle about points p1,p2,p3
c				using the algorithm of Cline and Renka
c				(see Sloan 1987).
c
c------------------------------------------------------------------------
c
	Function swap(p1,p2,p3,p)

	logical		swap

	real*8		p(2),p1(2),p2(2),p3(2)
	real*8		x13,y13,x23,y23,x1p,y1p,x2p,y2p
	real*8		cosa,cosb,sina,sinb

	x13=p1(1)-p3(1)
	y13=p1(2)-p3(2)
	x23=p2(1)-p3(1)
	y23=p2(2)-p3(2)
	x1p=p1(1)-p(1)
	y1p=p1(2)-p(2)
	x2p=p2(1)-p(1)
	y2p=p2(2)-p(2)

	cosa = x13*x23 + y13*y23
	cosb = x2p*x1p + y1p*y2p

	if((cosa.ge.0.d0).and.(cosb.ge.0.d0))then
            swap = .false.
        else if((cosa.lt.0.d0).and.(cosb.lt.0.d0))then
            swap = .true.
        else
            sina=x13*y23-x23*y13
            sinb=x2p*y1p-x1p*y2p
	    if((sina*cosb+sinb*cosa).lt.0.d0)then
                swap = .true.
            else
                swap = .false.
            end if
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c	Triloc_del - locates the triangle containing point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c		points(2,np)		array of node co-ordinates	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c	        loc			first guess of triangle containing
c					(x, y).
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        out			=true if (x,y) is outside of
c					the convex hull, otherwise = false. 
c	        k			index of face through which the
c					algorithm last passed (used by
c					routine visbilelist if out = .true.)
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a visible triangle
c		 on the hull, out is set to .true., and k is set to the
c		 index of the face of triangle loc visible from the input point
c		 (used as a starting point by the routine visiblelist)
c
c		 This version also returns the parameter iface. 
c		 If iface .ne. 0 then the input point is on the face of 
c		 triangle t between nodes iface and mod(iface,3)+1 and 
c		 it is also on the convex hull.
c
c		 A point is assumed to be on the edge (or its extension)
c		 between two nodes if it is inside the triangle at a 
c		 distance >= eps.
c
c		 Can be extended to higher dimensions using a similar
c		 stepping mechanism but without angular test.
c
c	         This routine was converted from using Sloan's definition of
c	         the adjacency matrix to the `opposite' definition on 30/1/96.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine Triloc_del
     &                       (x,y,points,vertices,neighbour,loc,eps,
     &                        out,k,iface)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        logical		out
	real*8		x,y,del1,del2,del
	real*8		eps
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/
        logical         new

	out = .false.
        new = .true.
        ic = 0

 10     continue
c					point is outside convex hull
        if(out)return
        iface = 0

        do 20 i=1,3
	   j = c1(i)
c	   k = c1(j)
c					definition of adjacency matrix

c					use Sloan's 
c					definition of adjacency matrix
	   k = i

           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
           del = del1-del2
 	   if(dabs(del).le.eps)then
              iface = i
	   else if(del.gt.0.d0)then
	      if(neighbour(c2(k),loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(c2(k),loc)
	      end if
              if(.not.new.and.loc.eq.loc1)then
                 write(*,100) 
                 write(*,*)' Current triangle:',loc, 
     &           ' last three:',loc1,loc2,loc3
                 write(*,*)' New point      x:',x,' y:',y
                 write(*,*)' Triangle ',loc,
     &           ' v:',(vertices(j,loc),j=1,3),
     &           ' n:',(neighbour(c2(j),loc),j=1,3)
c                write(*,*)' del',del,' del1',del1,' del2',del2
                 write(*,101) 
                 stop
              end if
              if(new)then
                ic = ic + 1
                if(ic.eq.3)new = .false.
              end if
              loc1 = loc2
              loc2 = loc3
              loc3 = loc
	      go to 10
	   end if
 20     continue
	
c						check if input point is
c						on the convex hull
c
        if(neighbour(c2(iface),loc).ne.0)iface = 0

c       if(iface.ne.0)then
c          j = mod(iface,3)+1
c          jj = vertices(iface,loc)
c          kk = vertices(j,loc)
c          write(*,*)' point on triangle between nodes ',
c    &               jj,' and',kk
c          write(*,*)' point is on the convex hull'
c       end if

 100    format(/'Error in subroutine Triloc_del:',//
     &  ' Infinite loop detected in walking triangle algorithm',/,
     &  ' Probably due to rounding error creating a flat triangle'/)

 101    format(/1x,'Remedy: '/
     &  ' Either increase size of parameter eps in calling routine '/
     &  ' or re-order input points by running program nn_hull '/)

	return
	end
c
c------------------------------------------------------------------------
c
c	Function ccw - used to test the orientation of three points
c
c		Input : 	points p1,p2 and p3 (vectors 2x1)
c				(e.g. p(1,2) = x co-ordinate of p2)
c
c		Output: 	ccw,I
c
c		ccw    k
c	  	 1     0	:The direction p1,p2,p3 is ccw (+ve)   
c	  	-1     0	:The direction p1,p2,p3 is  cw (-ve)   
c	  	 1     1	:p1,p2,p3 are colinear & p2 in middle  
c	  	-1     1	:p1,p2,p3 are colinear & p1 in middle
c	  	 0     1	:p1,p2,p3 are colinear & p3 in middle 
c
c
c				Calls no other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
      Integer Function ccw(p1,p2,p3,k)
c     
      real*8		p1(2),p2(2),p3(2)
      real*8		dx1,dx2,dy1,dy2,a,b
c     integer		ccw
c
      dx1 = p2(1) - p1(1)
      dx2 = p3(1) - p1(1)
      dy1 = p2(2) - p1(2)
      dy2 = p3(2) - p1(2)
      a = dx1*dy2
      b = dy1*dx2
      if (a.gt.b)then
         k = 0
         ccw = 1
      else if(a.lt.b)then
         k = 0
         ccw = -1
      else if(dx1*dx2.lt.0.0.or.dy1*dy2.lt.0.0)then
         k = 1
         ccw = -1 
      else if((dx1*dx1+dy1*dy1).lt.(dx2*dx2+dy2*dy2))then
         k = 1
         ccw = 1 
      else
         k = 1
         ccw = 0 
      end if
      return
      end
c
c------------------------------------------------------------------------
c
c	function theta - returns a real number between 0 and 360
c		         which has the same ordering as the angle
c			 between the line (a,b) and the horizontal.
c
c------------------------------------------------------------------------
c
	function theta(a,b)

	real*8		a(2),b(2)
	real*4		theta

	dx = b(1) - a(1)  
	ax = abs(dx)
	dy = b(2) - a(2)  
	ay = abs(dy)

	theta = 0.
	d = ax+ay
	if(d.ne.0)theta=dy/d

	if(dx.lt.0.0)then
           theta = 2.-theta
	else if(dy.lt.0.0)then
           theta = 4.+theta
	end if
	theta = theta*90.

	return
	end
c
c-----------------------------------------------------------------------------
c
c	Subroutine qhullf (np,i,j,nt_max,k,points,nt,vertices)
c
c-----------------------------------------------------------------------------
c
	Subroutine qhullf (np,i,j,nt_max,k,points,nt,vertices)
c
	real*8		points(2,*)
	integer		vertices(3,*)

	write(*,*)' '
	write(*,*)' Error in subroutine nn2d_setup'
	write(*,*)' qhull is not installed'
	write(*,*)' Delaunay triangulation must be either'
	write(*,*)' calculated with routine delaun (dmode=-1 or -2)'
	write(*,*)' or read in from a file (dmode>0; logical unit=dmode)'
	write(*,*)' '

	stop
	end
c
c------------------------------------------------------------------------
c
c       plot_c - dummy routine
c
c
c------------------------------------------------------------------------
c
        Subroutine plot_c(xs,ys,xs2,ys2)

	real*8		xs2,ys2
c                                               do nothing
        return
        end

c
c------------------------------------------------------------------------
c
c       plot_tc - dummy routine
c
c------------------------------------------------------------------------
c
        Subroutine plot_tc(n,points,vertices,centres)

        real*8          points(2,*)
        real*8          centres(3,*)
        integer         vertices(3,*)

c                                               do nothing
        return
        end
 
c
c ----------------------------------------------------------------------------
c
c       cputime - calls system dependent routine to calculate cputime
c		  since last call.
c
c       Calls dtime.
c						M. Sambridge, June 2001
c
c ----------------------------------------------------------------------------
c
        Function cputime(t1,t2)
        real*4 t1,t2
        real*4 tarray(2)
        
c        cputime = dtime(tarray)
c        t1 = tarray(1)
c        t2 = tarray(2)
        cputime=1.0

        return
        end
c------------------------------------------------------------------------
c
c	Ccentres - calculates centres of all Delaunay circumcircles
c
c
c	Input:
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c
c	Output:
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c
c	Comments:
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine ccentres(points,vertices,nt,centres)
c
	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x1,x2,x3,y1,y2,y3,x,y
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	integer		vertices(3,*)
c						Find centres of all
c						Delaunay Circumcircles
	do 5 i= 1,nt

	   x1 = points(1,vertices(1,i))
	   x2 = points(1,vertices(2,i))
	   x3 = points(1,vertices(3,i))
	   y1 = points(2,vertices(1,i))
	   y2 = points(2,vertices(2,i))
	   y3 = points(2,vertices(3,i))

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1
	   x = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &         (denom)

	   y = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/ 
     &         (denom)

	   centres(1,i) = x
	   centres(2,i) = y
           x1 = x - x1
           y1 = y - y1
	   centres(3,i) = x1*x1 + y1*y1

 5	continue

	return
	end
c
c------------------------------------------------------------------------
c
c	build_nv - Builds neighbour array for Delaunay triangulation in 2-D.
c
c	Input:	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c		np_max			maximum number of nodes
c		nmax			maximum total number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c
c	Output:
c		neighbour(3,nt)		array of neighbouring triangles
c
c	Comments:
c		 Assumes input list of vertices in anticlockwise sequence
c		 and produces an anticlockwise list of neighbour triangles.
c		 The value of neighbour(i,j) is the index of the neighbouring
c		 triangle opposite node i in triangle j.
c
c		 Three temporary work arrays are used and must be dimensioned
c		 in the calling program in the following way:
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntwork(nmax) : dummy array (NOTE NOT triangles
c					  attached to each node)
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c					(using ideas by J.Braun)
c
c------------------------------------------------------------------------
c
	Subroutine build_nv
     &             (np,vertices,nt,np_max,nmax,
     &              neighbour,nnn,nnlist,ntwork)
c
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntwork(*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud

	if(nnwrite)write(*,*)' Building neighbour v ...'
c
c					initialize neighbour list
	do 5 i = 1,3
	   do 4 j = 1,nt
	      neighbour(i,j) = 0
 4         continue
 5      continue
c					initialize work arrays
        do 6 i = 1,nmax
	   nnlist(i) = 0
	   ntwork(i) = 0
 6      continue

        do 7 i = 1,np
	   nnn(i) = 0
 7      continue

	do 10 it = 1,nt
	   i1 = vertices(1,it)
	   i2 = vertices(2,it)
	   i3 = vertices(3,it)
	   nnn(i1) = nnn(i1) + 1
	   nnn(i2) = nnn(i2) + 1
	   nnn(i3) = nnn(i3) + 1
 10     continue

c					turn nnn into a running sum
	itemp = nnn(1)+1
	nnn(1) = 1
	do 20 j = 2,np+1
	   itemp2  = itemp 
	   itemp   = itemp + nnn(j)+1
	   nnn(j) = itemp2 + 1
 20     continue
c       write(*,*)' size of array =',nnn(np+1)-1
c       write(*,*)' 3nt+np        =',3*nt+np

	if(nnn(np+1).ge.nmax)then
           write(*,*)'Error: array sizes too small in subroutine '
     &               ,'build_nv'
           write(*,*)'       maximum number of neighbours for all nodes'
           write(*,*)'       is too small: current value =',nmax
           write(*,*)'       Increase size of parameter nmax'
           write(*,*)'       to at least',nnn(np+1)
           write(*,*)'       This will be satisfied if nmax is set'
           write(*,*)'       to 3*nt_max+np_max in calling program' 
	   stop
	end if

	do 25 it = 1,nt
	   i1 = vertices(1,it) 
	   i2 = vertices(2,it) 
	   i3 = vertices(3,it) 
c						compare neighbours i1 i2
c						(remove go to ?)
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   jt = 0
	   do 30 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i2
	         ntwork(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 31
	      else if(nnlist(j).eq.i2.and.ntwork(j).ne.it)then
                 jt = ntwork(j)
	         go to 31
	      end if
  30       continue
  31       continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i2)
	      j2 = nnn(i2+1) - 1 
	      do 32 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntwork(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 33
	         end if
  32          continue
	   end if
  33       continue

	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i2
	      neighbour(3,it) = jt
	      k1 = vertices(1,jt)
	      k2 = vertices(2,jt)
	      k3 = vertices(3,jt)
	      if(k1.ne.i1.and.k1.ne.i2)then 
	         neighbour(1,jt) = it
	      else if(k2.ne.i1.and.k2.ne.i2)then 
	         neighbour(2,jt) = it
	      else
	         neighbour(3,jt) = it
	      end if
           end if
c						compare neighbours i1 i3
	   jt = 0
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   do 130 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntwork(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 131
	      else if(nnlist(j).eq.i3.and.ntwork(j).ne.it)then
                 jt = ntwork(j)
	         go to 131
	      end if
  130      continue
  131      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 132 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntwork(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 133
	         end if
  132         continue
	      end if
  133      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i3
	     neighbour(2,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i1.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i1.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if
c						compare neighbours i2 i3
	   jt = 0
	   j1 = nnn(i2)
	   j2 = nnn(i2+1) - 1 
	   do 230 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntwork(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 231
 	      else if(nnlist(j).eq.i3.and.ntwork(j).ne.it)then
                 jt = ntwork(j)
 	         go to 231
	      end if
  230      continue
  231      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 232 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i2
	            ntwork(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 233
	         end if
  232         continue
	      end if
  233      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i2 and i3
	     neighbour(1,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i2.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i2.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if

 25     continue

	if(nnwrite)write(*,*)' built neighbour v'

	return
	end
c------------------------------------------------------------------------
c
c       Calculate_hulltriangles - finds all triangles with a face on
c                                 the convex hull by searching through
c                                 the entries in the array neighbour.
c
c       Input:
c               neighbour(3,nt)         array of neighbouring tetrahedra
c               nt                      number of tetrahedra
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		nohalt			determines error response
c       Output:
c		hulltriangles(nh)	array of triangles with an edge
c					on the convex hull
c               nh                      number of tetrahedra with an edge
c                                       on the convex hull
c       Comments:
c
c                This routine fills up the array hulltriangles which
c                is only used by routine nn2Do, i.e the `pseudo-extension' 
c		 Watson's nn-interpolation method to points outside of the 
c		 convex hull. If nnext is set to false then hulltriangles
c		 is never used and the array can be set to size 1.
c
c		 If nohalt = 0 then the routine will stop with an error
c		 message if nh > nh_max. If nohalt .ne. 0 and nh > nh_max
c		 then it will return nh = -1. 
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, May 1995.
c
c------------------------------------------------------------------------
c
	Subroutine calculate_hulltriangles
     &             (neighbour,nt,nh_max,hulltriangles,nh,nohalt)
c
	integer		neighbour(3,*)
	integer		hulltriangles(*)

c                                               store list of triangles
c                                               which have an edge on the
c                                               convex hull.
c                                               (used by routine nn2D)
        nh = 1
        do 100 j = 1,nt
           if(neighbour(1,j).eq.0.or.
     &        neighbour(2,j).eq.0.or.
     &        neighbour(3,j).eq.0)then
              hulltriangles(nh) = j
              nh = nh + 1
              if(nh.gt.nh_max.and.nohalt.eq.0)then
                  write(*,*)' Error array storing outward facing '
                  write(*,*)' triangles on convex hull is too small.'
                  write(*,*)' Increase size of parameter nh_max'
                  stop
              else if(nh.gt.nh_max.and.nohalt.ne.0)then
                  nh = -1
                  return
              end if
           end if
 100    continue
        nh = nh -1
 
	return
	end
c
      program setbrn
c
c Herewith the new version of setbran.f with separated table and header
c and correct index assignments for direct access (hopefully)
c
      include 'limits.inc'
      save
      character*8 code,phcd
      character*20 modnam
      double precision zm,pm,pb,pu,taup,xp,taul,px,xt,xl,pux,pt,taut,
     1 coef,xa
      double precision tmp(nsl1,2),xm(nsl1,2),deg,dtol,zmax,zoc,zic,
     1 z0
      dimension ndx2(nsr0,2)
      common/umodc/zm(nsr0,2),pm(nsr0,2),ndex(nsr0,2),mt(2)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/brnc/phcd(jbrn)
      common/msqc/pux(jbrn,2),km(2),midx(jbrn,2)
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data nin,nout,xmin,dtol/1,2,200.,1d-6/
      deg=180d0/3.1415927d0
c
c     write(6,*) "rec length for dasign:"
c     read(5,*) ndasr
c
      call assign(nin,-1,'remodl.hed')
      read(nin)ndasr,modnam,zmax,zoc,zic,kb,(pb(i),i=1,kb(2)),
     1 mt,lt,lbb,lcb,xn,pn,tn
      read(nin)((lbrk(i,nph),i=1,lbb(nph)),(code(i,nph),
     1 i=1,lcb(nph)),(zm(i,nph),pm(i,nph),ndex(i,nph),i=1,mt(nph)),
     2 (lvz(i,nph),taul(i,nph),xl(i,nph),i=1,lt(nph)),nph=1,2)
      call retrns(nin)
      print *,'ndasr =',ndasr,'  modnam = ',modnam
c
      call dasign(nin,-1,'remodl.tbl',ndasr)
      nrec=0
      do 1 nph=1,2
      n1=kb(nph)
      ind=0
      do 2 k=1,n1
 2    xm(k,nph)=0d0
 3    nrec=nrec+1
      read(nin,rec=nrec)z0,n,(tmp(k,1),k=1,n),(tmp(k,2),k=1,n)
      if(ind.gt.0) go to 4
      if(dabs(z0-zoc).le.dtol) go to 4
      j=1
      do 5 k=2,n
      xm(k,nph)=dmax1(xm(k,nph),dabs(tmp(j,2)-tmp(k,2)))
 5    j=k
      if(n+1.eq.n1) xm(n1,nph)=tmp(n,2)
      go to 3
 4    ind=ind+1
      do 6 k=1,n
      taup(k,ind,nph)=tmp(k,1)
 6    xp(k,ind,nph)=tmp(k,2)
      if(ind.lt.3) go to 3
 1    continue
c
      xmin=xn*xmin
c
c^      call assign(10,2,'setbrn1.lis')
c^      write(10,*)'kb mt lt lbb lcb',kb,mt,lt,lbb,lcb
c^      write(10,*)'xn pn tn xmin',xn,pn,tn,xmin
      cn=1./xn
c^      write(10,209)(i,(lbrk(i,j),code(i,j),j=1,2),i=1,lbb(1))
 209  format(/(1x,2i5,2x,a,i5,2x,a))
c^      write(10,210)(i,lbrk(i,2),code(i,2),i=lbb(1)+1,lbb(2))
 210  format(1x,i5,15x,i5,2x,a)
c^      write(10,200,iostat=ios)(i,(zm(i,j),pm(i,j),ndex(i,j),j=1,2),
c^     1 i=1,mt(1))
 200  format(/(1x,i5,2f12.6,i5,2x,2f12.6,i5))
c^      write(10,201,iostat=ios)(i,zm(i,2),pm(i,2),ndex(i,2),
c^     1 i=mt(1)+1,mt(2))
 201  format(1x,i5,31x,2f12.6,i5)
c^      write(10,217)((nph,i,lvz(i,nph),taul(i,nph),deg*xl(i,nph),
c^     1 i=1,lt(nph)),nph=1,2)
 217  format(/(1x,3i5,f12.6,f12.2))
c^      write(10,202)(i,pb(i),cn*xm(i,1),cn*xm(i,2),i=1,kb(1))
 202  format(/(5x,i5,f12.6,2f12.2))
c^      write(10,203)(i,pb(i),cn*xm(i,2),i=kb(1)+1,kb(2))
 203  format(5x,i5,f12.6,12x,f12.2)
c^      call retrns(10)
c^      call assign(10,2,'setbrn2.lis')
c
      do 8 nph=1,2
      n1=kb(nph)
      do 9 i=2,n1
      xm(i,nph)=xm(i-1,nph)+xm(i,nph)
      pu(i,nph)=pb(i)
 9    kuse(i,nph)=-1
      do 8 ind=3,2,-1
      jnd=ind-1
      do 8 i=1,n1
      taup(i,ind,nph)=taup(i,ind,nph)-taup(i,jnd,nph)
 8    xp(i,ind,nph)=xp(i,ind,nph)-xp(i,jnd,nph)
      do 10 nph=1,2
 10   call pdecx(kb(nph),nph,xm,2.)
c
c^      write(10,*)'ku',ku
c^      write(10,204)(i,(pu(i,nph),cn*xm(i,nph),cn*(xm(i+1,nph)-
c^     1 xm(i,nph)),nph=1,2),i=1,ku(1))
 204  format(/(5x,i5,2(f12.6,2f12.2)))
c^      write(10,205)(i,pu(i,2),cn*xm(i,2),cn*(xm(i+1,2)-xm(i,2)),
c^     1 i=ku(1)+1,ku(2))
 205  format(5x,i5,36x,f12.6,2f12.2)
c^      do 207 nph=1,2
c^ 207  write(10,206)(i,pb(i),(taup(i,j,nph),j=1,3),(deg*xp(i,j,nph),
c^     1 j=1,3),i=1,kb(nph))
c^ 206  format(/(1x,i5,4f10.6,3f10.2))
c
      call layout
c     write(10,214)(i,pb(i),(kuse(i,j),j=1,2),i=1,kb(2))
c214  format(/(5x,i5,f12.6,2i5))
      do 11 nph=1,2
      n1=kb(nph)
      k=0
      do 12 i=1,n1
      if(kuse(i,nph).lt.0) go to 12
      k=k+1
      pu(k,nph)=pb(i)
 12   continue
 11   ku(nph)=k
      call kseq
      call mseq
c^      write(10,215)(i,(pu(i,j),j=1,2),i=1,ku(1))
 215  format(/(5x,i5,2f12.6))
c^      write(10,216)(i,pu(i,2),i=ku(1)+1,ku(2))
 216  format(5x,i5,12x,f12.6)
c^      write(10,208)(i,(nafl(i,j),j=1,3),(indx(i,j),j=1,2),(kndx(i,j),
c^     1 j=1,2),(fcs(i,j),j=1,3),i=1,nseg)
 208  format(/(1x,8i6,3f6.1))
c^      write(10,211)(i,(jndx(i,j),j=1,2),(mndx(i,j),j=1,2),(px(i,j),
c^     1 j=1,2),(deg*xt(i,j),j=1,2),phcd(i),i=1,nbrn)
 211  format(/(1x,i3,4i5,2f12.6,2f10.2,2x,a))
c^      write(10,218)(i,(midx(i,j),j=1,2),(pux(i,j),j=1,2),
c^     1 i=1,max0(km(1),km(2)))
 218  format(/(1x,i3,2i5,2f12.6))
c^      write(10,212,iostat=ios)(i,pt(i),taut(i),deg*xa(i),
c^     1 cn*(xa(i)-xa(i+1)),(coef(j,i),j=1,5),i=1,nl)
 212  format(/(1x,i4,0p2f12.6,2f10.2,1p5d10.2))
c^      call retrns(10)
c
c^      call assign(10,2,'setbrn3.lis')
      do 20 nph=1,2
      mt(nph)=mt(nph)-3
      ku(nph)=ku(nph)-1
 20   km(nph)=km(nph)-1
c     icor=33  -  originally 32 records used as header in setbrn
c                 and 2 records used as header in remodl.
      icor=3
      do 14 nph=1,2
      m1=mt(nph)
      icor=icor-3
      do 14 i=2,m1
 14   ndx2(i,nph)=ndex(i,nph)+icor
      len1=ku(2)+km(2)
      len0=8*len1
      len2=5*nl
c^      write(10,*)'nseg nbrn mt ku km len len1',nseg,nbrn,mt,ku,km,len0,
     1 len1
c^      write(10,*)
      nasgr = len0
      write(6,*) 'reclength for direct access', nasgr
c++ 
c     write(6,*) 'enter model name'
c     read(5,*) cmodel
      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)
c     cnam1 = cmodel(1:nb)//'.tbl'
c     cnam2 = cmodel(1:nb)//'.hed'
      write(6,*) 'header file  :',modnam(1:nb)//'.hed'
      write(6,*) 'table file   :',modnam(1:nb)//'.tbl'
      call assign(nout,-2,modnam(1:nb)//'.hed')
c++ 
      write(nout) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km,fcs,nafl,
     1 indx,kndx
      write(nout) pm,zm,ndx2
      write(nout) pu,pux
      write(nout) phcd,px,xt,jndx
      write(nout) pt,taut
      write(nout) coef
      call retrns(nout)
c
      call dasign(nout,-2,modnam(1:nb)//'.tbl',nasgr)
      nrec = 0
      do 16 nph=1,2
      m1=mt(nph)
      n1=ku(nph)
      k1=km(nph)
c^      write(10,*)'nph m1 n1 k1',nph,m1,n1,k1
      do 16 m=2,m1
      if(ndex(m,nph).eq.ndex(m-1,nph)) go to 16
      read(nin,rec=ndex(m,nph))z0,n,(tmp(k,1),k=1,n),(tmp(k,2),k=1,n)
c^      write(10,*)'m nph ndex n',m,nph,ndex(m,nph),n
      k=0
      l=1
      do 17 i=1,n
      if(kuse(i,nph).lt.0) go to 17
      if(dabs(pux(l,nph)-pb(i)).gt.dtol) go to 18
      tmp(l,2)=tmp(i,2)
      l=l+1
 18   k=k+1
      tmp(k,1)=tmp(i,1)
 17   continue
c^      write(10,*)'k l nrec',k,l-1,nrec+1,ndx2(m,nph),sngl(tmp(1,1))
      if(k.ge.n1) go to 19
      k=k+1
      do 21 i=k,n1
 21   tmp(i,1)=0d0
 19   if(l.gt.k1) go to 23
      do 22 i=l,k1
 22   tmp(i,2)=0d0
 23   nrec=nrec+1
      write(nout,rec=nrec)(tmp(i,1),i=1,n1),(tmp(i,2),i=1,k1)
 16   continue
c
      call retrns(10)
      call retrns(nin)
      call retrns(nout)
      call vexit(0)
      end
c
      subroutine pdecx(n1,nph,xm,fac)
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl
      double precision xm(nsl1,2),ptol,pa,pax,plim
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      data ptol/.03/
c
      call collct(1,n1,xm(1,nph),fac*xmin)
      k=0
      plim=.7d0*pu(n1,nph)
      do 1 i=1,n1
      if(xm(i,nph).lt.0d0) go to 1
      if(pu(i,nph).lt.plim) go to 2
      if(pu(i,nph)-pu(k,nph).le.ptol) go to 2
      pa=pu(k,nph)+.75d0*(pu(i,nph)-pu(k,nph))
      pax=1d10
      m=0
      do 3 j=i1,i
      if(dabs(pu(j,nph)-pa).ge.pax) go to 3
      m=j
      pax=dabs(pu(j,nph)-pa)
 3    continue
      if(m.eq.i1.or.m.eq.i) go to 2
      k=k+1
      pu(k,nph)=pu(m,nph)
      xm(k,nph)=0d0
      kuse(m,nph)=1
 2    k=k+1
      i1=i
      pu(k,nph)=pu(i,nph)
      xm(k,nph)=xm(i,nph)
      kuse(i,nph)=1
 1    continue
      ku(nph)=k
      return
      end
c
      subroutine collct(i1,i2,x,xmn)
c
c $$$$$ calls varn $$$$$
c
      save
      double precision x(i2)
      data cn/6371./
c
      is=i1+1
      ie=i2-1
      if(ie.lt.is) return
      k1=i1
      var=0.
      m=0
      do 1 i=is,ie
      dx1=dabs(x(k1)-x(i))-xmn
      dx2=dabs(x(k1)-x(i+1))-xmn
      if(abs(dx2).ge.abs(dx1)) go to 2
      x(i)=-x(i)
      go to 1
 2    if(k1.le.i1) kb=i
      k1=i
      var=var+dx1*dx1
      m=m+1
 1    continue
      dx1=dabs(x(k1)-x(i2))-xmn
      var=var+dx1*dx1
      m=m+1
 7    if(m.le.1) return
      k1=i1
      k2=kb
      ks=kb+1
      nch=0
      do 8 i=ks,i2
      if(x(i).lt.0d0) go to 8
      k0=k1
      k1=k2
      k2=i
      var1=varn(x,k0,k1,k2,k1-1,xmn,var,m,m1)
      var2=varn(x,k0,k1,k2,k1+1,xmn,var,m,m2)
      if(amin1(var1/m1,var2/m2).ge.var/m) go to 6
      nch=nch+1
      x(k1)=-x(k1)
      if(var1/m1-var2/m2)3,4,5
 4    if(m1-m2)3,3,5
 3    k1=k1-1
      x(k1)=dabs(x(k1))
      var=var1
      m=m1
      go to 6
 5    k1=k1+1
      x(k1)=dabs(x(k1))
      var=var2
      m=m2
 6    if(k0.eq.i1) kb=k1
 8    continue
      if(nch.gt.0) go to 7
      return
      end
c
      function varn(x,k0,k1,k2,kt,xmn,var,m,mn)
c
c $$$$$ calls only library routines $$$$$
c
      save
      double precision x(k2)
c
      dx1=dabs(x(k0)-x(k1))-xmn
      dx2=dabs(x(k1)-x(k2))-xmn
      varn=var-dx1*dx1-dx2*dx2
      if(kt.le.k0.or.kt.ge.k2) go to 1
      dx1=dabs(x(k0)-dabs(x(kt)))-xmn
      dx2=dabs(dabs(x(kt))-x(k2))-xmn
      varn=varn+dx1*dx1+dx2*dx2
      mn=m
      return
 1    dx1=dabs(x(k0)-dabs(x(k2)))-xmn
      varn=varn+dx1*dx1
      mn=m-1
      return
      end
c
      subroutine layout
c
c   Layout contains the program for the desired travel-time segments 
c   implemented as calls to the mk_br entry points.  Each call does one 
c   segment (which may have many branches).
c
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision dir(3),cref(3),sref(3)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data dir,cref,sref/1d0,1d0,1d0,1d0,2d0,2d0,2d0,2d0,2d0/
c
c   Initialize variables.
      nseg=0
      nbrn=0
      nl=0
      do 1 j=1,3
      do 1 i=1,jseg
 1    fcs(i,j)=0.
      do 2 i=1,jout
 2    taut(i)=0d0
      do 3 j=1,2
      do 3 i=1,jbrn
 3    xt(i,j)=0d0
c
c   Do all of the segments.
c
c   P (up-going branch).
      print *,'Layout:  do Pup'
      call mkubr(ku(1),   +1)
c   P, Pdiff, and PKP.
      print *,'Layout:  do P and PKP'
      call mkdbr(1,lbb(1),-1,3,1,1,dir)
c   PKiKP.
      print *,'Layout:  do PKiKP'
      call mkrbr(2,       -1,2,1,1,dir)
c   pP.
      print *,'Layout:  do pP'
      call mkdbr(1,lbb(1),+1,3,1,1,dir)
c   sP.
      print *,'Layout:  do sP'
      call mkdbr(1,lbb(1),+2,3,1,1,dir)
c   pPKiKP.
      print *,'Layout:  do pPKiKP'
      call mkrbr(2,       +1,2,1,1,dir)
c   sPKiKP.
      print *,'Layout:  do sPKiKP'
      call mkrbr(2,       +2,2,1,1,dir)
c   PcP.
      print *,'Layout:  do PcP'
      call mkrbr(3,       -1,1,1,1,dir)
c   ScP.
      print *,'Layout:  do ScP'
      call mkrbr(3,       -2,1,2,1,dir)
c   SKP.
      print *,'Layout:  do SKP'
      call mkdbr(1,3,     -2,3,2,1,dir)
c   SKiKP.
      print *,'Layout:  do SKiKP'
      call mkrbr(2,       -2,2,2,1,dir)
c   PKKP.
      print *,'Layout:  do PKKP'
      call mkdbr(1,3,     -1,3,1,1,cref)
c   SKKP.
      print *,'Layout:  do SKKP'
      call mkdbr(1,3,     -2,3,2,1,cref)
c   PP and P'P'.
      print *,'Layout:  do PP, P''P'''
      call mkdbr(1,lbb(1),-1,3,1,1,sref)
c   S (up-going branch).
      print *,'Layout:  do Sup'
      call mkubr(ku(2),   +2)
c   S, Sdiff, and SKS.
      print *,'Layout:  do S and SKS'
      call mkdbr(1,lbb(2),-2,3,2,2,dir)
c   pS
      print *,'Layout:  do pS'
      call mkdbr(1,lbb(1),+1,3,2,2,dir)
c   sS
      print *,'Layout:  do sS'
      call mkdbr(1,lbb(2),+2,3,2,2,dir)
c   ScS
      print *,'Layout:  do ScS'
      call mkrbr(4,       -2,1,2,2,dir)
c   PcS
      print *,'Layout:  do PcS'
      call mkrbr(3,       -1,1,1,2,dir)
c   PKS
      print *,'Layout:  do PKS'
      call mkdbr(1,3,     -1,3,1,2,dir)
c   PKKS
      print *,'Layout:  do PKKS'
      call mkdbr(1,3,     -1,3,1,2,cref)
c   SKKS
      print *,'Layout:  do SKKS'
      call mkdbr(1,3,     -2,3,2,2,cref)
c   SS and S'S'.
      print *,'Layout:  do SS and S''S'''
      call mkdbr(1,lbb(2),-2,3,2,2,sref)
c   SP
      print *,'Layout:  do SP'
      call mkcbr(4,lbb(1),-2,1,2,1,sref)
c   PS
      print *,'Layout:  do PS'
      call mkcbr(4,lbb(1),-1,1,1,2,sref)
      return
      end
c
      subroutine mkdbr(l1,l2,isgn,lyr,nph,kph,fac)
c
c   Mkdbr sets up a simple refracted wave segment.  L1 and l2 point to the 
c   lbrk array of slowness break point pointers.  Note that the P and S 
c   break point arrays don't necessarily line up layer by layer.  This is 
c   not generally a problem as most phases need only worry about the 
c   pointer to the surface slowness for one wave type and a pointer 
c   somewhere in the core (which is constrained to be the same for both 
c   P and S).  Isgn is positive if the wave starts out going up and 
c   negative if the wave starts out going down.  Iabs(isng) is 1 if the 
c   wave starts out as a P wave and 2 if the wave starts out as an S wave.  
c   Lyr gives the number of major layers (mantle, outer core, and inner 
c   core) that the wave penetrates.  Nph and kph give the wave type (1 for 
c   P and 2 for S) on the down-going and up-going legs of the ray path 
c   respectively.  Fac is a three element array giving the number of 
c   repeats of the ray path in each major layer.  This scheme incorporates 
c   turning rays (e.g., P and S), turning rays reflected, but not 
c   converted at the surface (e.g., PP and SS), up-going rays reflected 
c   and/or converted at the surface into turning rays (e.g., pP and sP), 
c   turning rays converted during transmission through an interface (e.g., 
c   SKP and PKS), and rays which turn multiple times while reflecting from 
c   the bottom side of a layer (e.g., PKKP or SKKP).  Mkdbr does not 
c   include up-going rays (to the receiver), rays reflected from the top 
c   side of a discontinuity, or rays which are reflected and converted at 
c   the free surface.  See mkubr, mkrbr, and mkcbr respectively for 
c   routines which handle these types of rays.
c
      save
      include 'limits.inc'
      character*8 code,phcd,ks
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision fac(3)
      common/brkptb/lcb(2),lbb(2),lbrk(nbr2,2)
      common/brkptc/code(nbr1,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/brnc/phcd(jbrn)
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
      data ks/'KKKKKKKK'/
c
c   Remember the programming as part of the final phase construction is 
c   done in depcor.
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
c   Using l1 and l2 to get the breakpoints has some shortcommings, 
c   particularly for converted phases.  It would be more general to 
c   have separate indicies for the breakpoints and the layers covered.
      if(l1.gt.1) kndx(nseg,1)=lbrk(l1-1,nph)
      kndx(nseg,2)=min0(min0(lbrk(l2,nph),lbrk(l2,kph)),
     1 lbrk(l2,iabs(isgn)))
      print *,'Mkdbr:  l1 l2 isgn lyr nph kph =',l1,l2,isgn,lyr,nph,kph
      print *,'Mkdbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 14 m=1,lyr
      fcs(nseg,m)=fac(m)
 14   xfc=amax1(xfc,fcs(nseg,m))
c
c   Set up the required slownesses, taus and distances.
c
      j=kndx(nseg,1)
      lz1=1
      lz2=1
c   Loop over the layers of interest.
      do 1 i=l1,l2
c   Be sure that the phase cuts off at the right place.
      l=min0(lbrk(i,nph),kndx(nseg,2))
c   Skip all total internal reflections.
      if(code(i,nph)(1:1).eq.'r'.or.j.ge.l) go to 1
c   Set the starting branch pointer.
      nbrn=nbrn+1
      nt=nl+1
      jndx(nbrn,1)=nt
c   Copy in the desired slownesses.
      do 2 k=j,l
      nl=nl+1
      pt(nl)=pb(k)
c   Add up the tau contributions.
      do 2 m=1,lyr
 2    taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
c   Take care of branch end pointers and slownesses.
      mndx(nbrn,1)=j
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(j)
      px(nbrn,2)=pb(l)
c   Add up distance contributions for the branch end points only.
      do 3 m=1,lyr
      xt(nbrn,1)=xt(nbrn,1)+fac(m)*(xp(j,m,nph)+xp(j,m,kph))
 3    xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
c   Take care of the contribution due to low velocity zones for the 
c   down-going leg(s).
      if(j.ne.lvz(lz1,nph)) go to 9
      do 11 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,nph)
 11   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,nph)
      taut(nt)=taut(nt)+fac(1)*taul(lz1,nph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz1,nph)
      lz1=lz1+1
c   Take care of the contributions due to low velocity zones for the 
c   up-going leg(s).
 9    if(j.ne.lvz(lz2,kph)) go to 10
      do 12 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,kph)
 12   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,kph)
      taut(nt)=taut(nt)+fac(1)*taul(lz2,kph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz2,kph)
      lz2=lz2+1
c   Decimate the slownesses if the branch is oversampled in distance.
 10   call pdect(jndx(nbrn,1),nl,j,iabs(isgn),xfc)
c   Set up the interpolation.
      call tauspl(jndx(nbrn,1),nl,pt,coef)
c   Remember the final branch end slowness value.
      jndx(nbrn,2)=nl
c
c   Take care of the branch name.  First, set up a default.
      phcd(nbrn)=code(i,nph)(2:2)//code(i,kph)(3:)
      if(idint(fac(1)+.5d0).gt.1) go to 5
      if(idint(fac(2)+.5d0).le.1) go to 4
c   Re-do the name if the ray is reflected from the underside of the 
c   core-mantle boundary.
      ind=idint(fac(2)-.5d0)
      phcd(nbrn)=code(i,nph)(2:2)//ks(1:ind)//code(i,kph)(3:)
      go to 4
c   Re-do the name if the ray is reflected from the surface.
 5    if(code(i,nph)(3:3).eq.' ') phcd(nbrn)=code(i,nph)(2:2)//
     1 code(i,kph)(2:)
      if(code(i,nph)(3:3).ne.' '.and.code(i,nph)(3:3).ne.'K')
     1  phcd(nbrn)=code(i,nph)(2:3)//code(i,kph)(2:)
      if(code(i,nph)(3:3).eq.'K') phcd(nbrn)=code(i,nph)(2:2)//''''//
     1 code(i,kph)(2:2)//''''//code(i,kph)(5:)
c   Take care .
 4    ind=max0(index(phcd(nbrn),'KSab'),index(phcd(nbrn),'S''ab'))
      if(phcd(nbrn)(1:1).eq.'S'.and.ind.gt.0) phcd(nbrn)(ind+2:ind+3)=
     1 'ac'
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
 1    j=l
      indx(nseg,2)=nl
      return
c
c   Mkubr handles up-going P and S.  L1 and isgn are as for mkdbr (except 
c   that l1 actually plays the role of l2 with the beginning break point 
c   assumed to be zero).  The other arguments are not needed.
c
      entry mkubr(l1,isgn)
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=0
      nafl(nseg,3)=0
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      l=kb(iabs(isgn))
      kndx(nseg,2)=l
      print *,'Mkubr:  l1 isgn =',l1,isgn
      print *,'Mkubr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      nbrn=nbrn+1
      jndx(nbrn,1)=nl+1
      do 6 k=1,l1
      nl=nl+1
      pt(nl)=pu(k,isgn)
 6    xa(nl)=0d0
      mndx(nbrn,1)=1
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(1)
      px(nbrn,2)=pb(l)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
      phcd(nbrn)=code(1,iabs(isgn))(2:2)
      indx(nseg,2)=nl
      return
c
c   Mkrbr handles reflected phases possibly with a conversion such as 
c   PcP, PcS, and PkiKP.  Arguments are as for mkdbr (except that l1 
c   actually plays the role of l2 with the beginning break point assumed 
c   to be zero).
c
      entry mkrbr(l1,isgn,lyr,nph,kph,fac)
      nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      l=min0(lbrk(l1,nph),lbrk(l1,kph))
      kndx(nseg,2)=l
      print *,'Mkrbr:  l1 isgn lyr nph kph =',l1,isgn,lyr,nph,kph
      print *,'Mkrbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 15 m=1,lyr
      fcs(nseg,m)=fac(m)
 15   xfc=amax1(xfc,fcs(nseg,m))
      if(lyr.ge.2) xfc=2.
c
      nbrn=nbrn+1
      jndx(nbrn,1)=nl+1
      do 7 k=1,l
      nl=nl+1
      pt(nl)=pb(k)
      do 7 m=1,lyr
 7    taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
      mndx(nbrn,1)=1
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(1)
      px(nbrn,2)=pb(l)
      do 8 m=1,lyr
 8    xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
      call pdect(jndx(nbrn,1),nl,1,iabs(isgn),xfc)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
      if(lyr.eq.1) phcd(nbrn)=code(l1,nph)(2:2)//'c'//code(l1,kph)(2:2)
      if(lyr.eq.2) phcd(nbrn)=code(l1,nph)(2:2)//code(l1,kph)(3:)
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
      indx(nseg,2)=nl
      return
c
c   Mkcbr handles phases reflected and converted at the surface such as 
c   PS and SP.  Arguments are as for mkdbr.
c
      entry mkcbr(l1,l2,isgn,lyr,nph,kph,fac)
      if(nph.gt.0.and.kph.gt.0.and.nph.ne.kph) go to 29
      print *,'Mkcbr:  bad call - nph kph =',nph,kph
      call vexit(1)
 29   nseg=nseg+1
      nafl(nseg,1)=isgn
      nafl(nseg,2)=nph
      nafl(nseg,3)=kph
      indx(nseg,1)=nl+1
      kndx(nseg,1)=1
      if(l1.gt.1) kndx(nseg,1)=min0(lbrk(l1,nph),lbrk(l1,kph))
      kndx(nseg,2)=min0(min0(lbrk(l2,nph),lbrk(l2,kph)),
     1 lbrk(l2,iabs(isgn)))
      print *,'Mkcbr:  l1 l2 isgn lyr nph kph =',l1,l2,isgn,lyr,nph,kph
      print *,'Mkcbr:  nseg kndx indx =',nseg,kndx(nseg,1),kndx(nseg,2),
     1 indx(nseg,1)
      xfc=0.
      do 16 m=1,lyr
      fcs(nseg,m)=fac(m)
 16   xfc=amax1(xfc,fcs(nseg,m))
c
      j=kndx(nseg,1)
      lz1=1
      lz2=1
      ik=l1
c
      print *,'Mkcbr:  start loop'
      do 17 in=l1,l2
 31   l=min0(lbrk(in,nph),kndx(nseg,2))
      if(code(in,nph)(1:1).eq.'r'.or.j.ge.l) go to 17
      l=min0(lbrk(ik,kph),kndx(nseg,2))
      if(code(ik,kph)(1:1).ne.'r'.and.j.lt.l.or.ik.ge.l2) go to 28
      j=max0(j,l)
      ik=ik+1
      go to 31
c
 28   if(lbrk(in,nph).le.lbrk(ik,kph)) go to 26
      l=min0(lbrk(ik,kph),kndx(nseg,2))
      print *,'kph:  kph ik j l code =',kph,ik,j,l,' ',code(ik,kph)
      isw=2
      go to 27
 26   l=min0(lbrk(in,nph),kndx(nseg,2))
      print *,'nph:  nph in j l code =',nph,in,j,l,' ',code(in,nph)
      isw=1
c
 27   nbrn=nbrn+1
      nt=nl+1
      jndx(nbrn,1)=nt
      do 18 k=j,l
      nl=nl+1
      pt(nl)=pb(k)
      do 18 m=1,lyr
 18   taut(nl)=taut(nl)+fac(m)*(taup(k,m,nph)+taup(k,m,kph))
      mndx(nbrn,1)=j
      mndx(nbrn,2)=l
      px(nbrn,1)=pb(j)
      px(nbrn,2)=pb(l)
      do 19 m=1,lyr
      xt(nbrn,1)=xt(nbrn,1)+fac(m)*(xp(j,m,nph)+xp(j,m,kph))
 19   xt(nbrn,2)=xt(nbrn,2)+fac(m)*(xp(l,m,nph)+xp(l,m,kph))
      if(j.ne.lvz(lz1,nph)) go to 20
      do 21 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,nph)
 21   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,nph)
      taut(nt)=taut(nt)+fac(1)*taul(lz1,nph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz1,nph)
      lz1=lz1+1
 20   if(j.ne.lvz(lz2,kph)) go to 22
      do 23 m=1,lyr
      taut(nt)=taut(nt)-fac(m)*taup(j,m,kph)
 23   xt(nbrn,1)=xt(nbrn,1)-fac(m)*xp(j,m,kph)
      taut(nt)=taut(nt)+fac(1)*taul(lz2,kph)
      xt(nbrn,1)=xt(nbrn,1)+fac(1)*xl(lz2,kph)
      lz2=lz2+1
 22   call pdect(jndx(nbrn,1),nl,j,iabs(isgn),xfc)
      call tauspl(jndx(nbrn,1),nl,pt,coef)
      jndx(nbrn,2)=nl
c
      if(code(in,nph)(3:3).eq.' ') phcd(nbrn)=code(in,nph)(2:2)//
     1 code(ik,kph)(2:)
      if(code(in,nph)(3:3).ne.' '.and.code(in,nph)(3:3).ne.'K')
     1  phcd(nbrn)=code(in,nph)(2:3)//code(ik,kph)(2:)
      if(code(in,nph)(3:3).eq.'K') phcd(nbrn)=code(in,nph)(2:2)//''''//
     1 code(ik,kph)(2:2)//''''//code(ik,kph)(5:)
      if(isgn.eq.1) phcd(nbrn)='p'//phcd(nbrn)
      if(isgn.eq.2) phcd(nbrn)='s'//phcd(nbrn)
      print *,'phcd:  in ik phcd =',in,ik,' ',phcd(nbrn)
      if(isw.le.1) go to 17
      ik=ik+1
      j=max0(j,l)
      go to 31
 17   j=max0(j,l)
      indx(nseg,2)=nl
      return
      end
c
      subroutine pdect(i1,i2,j1,nph,fac)
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pt,taut,coef,xa
      double precision h1,h2,hh
      dimension ib(2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/outc/pt(jout),taut(jout),coef(5,jout),xa(jout),nl
c
      xmn=fac*xmin
      isg=1
      do 1 i=1,2
      ib(i,1)=i1
 1    ib(i,2)=i2
      ii=i1+1
      ie=i2-1
      xa(i1)=xt(nbrn,1)
      do 2 i=ii,ie
      h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
 2    xa(i)=-(h2*taut(i-1)-(h2+h1)*taut(i)+h1*taut(i+1))/hh
      xa(i2)=xt(nbrn,2)
      do 3 i=ii,ie
      if((xa(i+1)-xa(i))*(xa(i)-xa(i-1)).gt.0d0) go to 3
      isg=2
      ib(1,2)=i-2
      ib(2,1)=i+2
 3    continue
      do 4 it=1,isg
 4    call collct(ib(it,1),ib(it,2),xa,xmn)
      k=i1-1
      j=j1
      do 5 i=i1,i2
      if(xa(i).lt.0d0) go to 5
      k=k+1
      pt(k)=pt(i)
      taut(k)=taut(i)
      xa(k)=xa(i)
      kuse(j,nph)=1
 5    j=j+1
      if(k.eq.nl) return
      ii=k+1
      do 6 i=ii,nl
 6    taut(i)=0d0
      nl=k
      return
      end
c
      subroutine kseq
c
c   Kseq makes a correspondence between model slownesses in array pb and 
c   the subset of the same slownesses used for sampling tau which are 
c   stored in pu (separate sets for P and S).  The net result is to 
c   translate the kndx pointers to critical slowness values (bounding 
c   branches actually implemented) from pointing into pb to pointing 
c   into pu.
c   
c
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl
      dimension kl(2),kk(jseg,2,2)
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      data kl/0,0/
c
c   Compile a sorted list of unique kndx values in the first column of 
c   kk.
c
      do 1 i=1,nseg
      nph=iabs(nafl(i,1))
      k=kl(nph)
      do 2 j=1,2
      if(k.le.0) go to 3
      do 4 m=1,k
      if(kk(m,nph,1)-kndx(i,j))4,2,5
 4    continue
 3    k=k+1
      kk(k,nph,1)=kndx(i,j)
      go to 2
 5    do 6 l=k,m,-1
 6    kk(l+1,nph,1)=kk(l,nph,1)
      k=k+1
      kk(m,nph,1)=kndx(i,j)
 2    continue
 1    kl(nph)=k
c
c   Make the correspondence between pb and pu for each kndx and save it 
c   in the second column of kk.
c
      do 7 nph=1,2
      n1=ku(nph)
      k=1
      ki=kk(k,nph,1)
      do 8 i=1,n1
      if(pu(i,nph)-pb(ki))8,9,10
 9    kk(k,nph,2)=i
      if(k.ge.kl(nph)) go to 7
      k=k+1
      ki=kk(k,nph,1)
 8    continue
 10   write(*,100)ki,pb(ki),nph
 100  format(1x,'Kseq:  pb(',i3,') =',f7.4,' not found in pu(*,',i1,
     1 ').')
      call vexit(1)
 7    continue
c
c   Replace each kndx pb pointer with the corresponding pu pointer.
c
      do 11 i=1,nseg
      nph=iabs(nafl(i,1))
      k=kl(nph)
      do 11 j=1,2
      do 12 m=1,k
      if(kk(m,nph,1)-kndx(i,j))12,11,13
 12   continue
 13   write(*,101)kndx(i,j)
 101  format(1x,'Kseq:  kndx value',i4,' not translated.')
      call vexit(1)
 11   kndx(i,j)=kk(m,nph,2)
      return
      end
c
      subroutine mseq
c         partial reordering of tables
      save 
      include 'limits.inc'
      double precision pb,pu,taup,xp,taul,xl,px,xt,pux
      common/intc/pb(nsl1),pu(nsl1,2),taup(nsl1,3,2),xp(nsl1,3,2),
     1 taul(nlvz0,2),xl(nlvz0,2),xmin,kb(2),ku(2),lt(2),lvz(nlvz0,2),
     2 kuse(nsl1,2)
      common/segc/fcs(jseg,3),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),
     1 nseg
      common/brnb/px(jbrn,2),xt(jbrn,2),jndx(jbrn,2),mndx(jbrn,2),nbrn
      common/msqc/pux(jbrn,2),km(2),midx(jbrn,2)
      data km/0,0/
c
      is=1
      do 1 i=1,nbrn
 8    if(jndx(i,2).le.indx(is,2)) go to 7
      is=is+1
      go to 8
 7    nph=iabs(nafl(is,1))
      k=km(nph)
      do 2 j=1,2
      if(k.le.0) go to 3
      do 4 m=1,k
      if(midx(m,nph)-mndx(i,j))4,2,5
 4    continue
 3    k=k+1
      midx(k,nph)=mndx(i,j)
      pux(k,nph)=px(i,j)
      go to 2
 5    do 6 l=k,m,-1
      midx(l+1,nph)=midx(l,nph)
 6    pux(l+1,nph)=pux(l,nph)
      k=k+1
      midx(m,nph)=mndx(i,j)
      pux(m,nph)=px(i,j)
 2    continue
 1    km(nph)=k
      return
      end
*-----        Location Routines      -------------------------------- + ----
c                                                                     ydist
      subroutine ydist
     ^        (dlats, dlons, dlatr, dlonr, delta, cazim, bazim, azima)
c
c     AUTHOR:  Brian L.N. Kennett  RSES, ANU
c     DATE:    January 1985
c     PURPOSE:
c             YDIST        Calculates distance and azimuth
c                          for spheroidal earth between
c                          specified geographic source and
c                          receiver station coordinates
c
*----------------------------------------------------------------------*------*
c     PARAMETERS
c
      real    dlats, dlons, dlatr, dlonr, delta, cazim, bazim
c
c     dlats  latitude of source
c     dlons  longitude of source
c     dlatr  latitude of receiver
c     dlonr  longitude of receiver
c     delta  angular distance
c     cazim  apparent azimuth at an array
c     bazim   azimuth from epicentre to receiver
c
*----------------------------------------------------------------------*------*
c
c     implicit real*8 (a-h,o-z)
      real*8 ecc,re,ec1,pi,pib2,degr,rlats,rlons,rlatr
      real*8 rlonr,glats,glatr,sps,cps,spr,cpr,rs,rr
      real*8 trs,prs,trr,prr,AS,BS,CS,DS,ES,GS,HS,KS
      real*8 AR,BR,CR,DR,ER,GR,HR,KR
      real*8 cosdr,deltar,sindr,deltak,szs,czs,szr,czr
      real*8 e,x,y
      real azima
c                          radius on spheroid
      gra(x,y,e) = dsqrt( (1.0d0-e)**2 /
     &                   ((1.0d0-e*y)**2 + e*e*x*y ) )
      ecc = 0.003367
      re = 6378.388
      ec1 = (1.0d0-ecc)**2
      pi = 3.141592653589793
      pib2 = pi/2.0
      degr = pi/180.0
      rlats = dlats*degr
      rlons = dlons*degr
      rlatr = dlatr*degr
      rlonr = dlonr*degr
c                          geocentric coordinates
      glats = datan2 ( ec1*dsin(rlats) ,dcos(rlats) )
      glatr = datan2 ( ec1*dsin(rlatr) ,dcos(rlatr) )
      sps = dsin(glats)**2
      cps = dcos(glats)**2
      spr = dsin(glatr)**2
      cpr = dcos(glatr)**2
c                          radii at source,receiver
      rs = re*gra(sps,cps,ecc)
      rr = re*gra(spr,cpr,ecc)
c
      trs = pib2 - glats
      prs = dlons*degr
      trr = pib2 - glatr
      prr = dlonr*degr
c                          direction cosines for source
      AS = dsin(trs)*dcos(prs)
      BS = dsin(trs)*dsin(prs)
      CS = dcos(trs)
      DS = dsin(prs)
      ES = -dcos(prs)
      GS = dcos(trs)*dcos(prs)
      HS = dcos(trs)*dsin(prs)
      KS = -dsin(trs)
c                          direction cosines for receiver
      AR = dsin(trr)*dcos(prr)
      BR = dsin(trr)*dsin(prr)
      CR = dcos(trr)
      DR = dsin(prr)
      ER = -dcos(prr)
      GR = dcos(trr)*dcos(prr)
      HR = dcos(trr)*dsin(prr)
      KR = -dsin(trr)
c                          distance
      cosdr = AS*AR + BS*BR + CS*CR
      deltar = dacos(cosdr)
      sindr = dsin(deltar)
c
      deltak = deltar*0.5d0*(rr+rs)
      delta = deltar/degr
c                          azimuth
      szs = DS*AR + ES*BR
      czs = GS*AR + HS*BR + KS*CR
      szr = DR*AS + ER*BS
      czr = GR*AS + HR*BS + KR*CS
c                          azima - azimuth to source
c                          bazim - backazimuth from source
c                          cazim - apparent azimuth at an array
      if (szr.eq.0.0) then
        bazim = 0.0
        if(dlats.gt.dlatr)then
           azima = 360.0
        else
           azima = 180.0
        endif
      else
        bazim = datan2(-szs ,-czs ) /degr
        azima = datan2(-szr ,-czr ) /degr
      end if
      if( bazim .lt. 0.0) bazim = bazim + 360.0
      cazim = azima + 180.0
      if( azima .lt. 0.0) azima = azima + 360.0
c
      if( cazim.lt. 0.0) cazim = cazim + 360.0
c
      return
      end
