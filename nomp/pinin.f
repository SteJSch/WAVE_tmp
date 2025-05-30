*CMZ :          06/05/2024  16.54.16  by  Michael Scheer
*CMZ :  4.01/04 22/11/2023  17.20.47  by  Michael Scheer
*CMZ :  4.00/15 12/02/2022  17.08.21  by  Michael Scheer
*CMZ :  4.00/13 04/12/2021  12.10.40  by  Michael Scheer
*CMZ :  4.00/01 05/04/2019  15.09.32  by  Michael Scheer
*CMZ :  3.06/00 25/02/2019  17.19.30  by  Michael Scheer
*CMZ :  3.02/08 24/06/2015  16.06.32  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.59.38  by  Michael Scheer
*CMZ :  2.67/02 28/03/2012  08.22.03  by  Michael Scheer
*CMZ :  2.66/06 22/05/2010  16.48.23  by  Michael Scheer
*CMZ :  2.64/06 15/09/2009  11.10.54  by  Michael Scheer
*CMZ :  2.64/05 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.62/03 16/07/2007  11.51.09  by  Michael Scheer
*CMZ :  2.62/02 16/07/2007  06.52.55  by  Michael Scheer
*CMZ :  2.54/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/10 05/11/2004  16.59.34  by  Michael Scheer
*CMZ :  2.51/00 17/05/2004  17.43.59  by  Michael Scheer
*CMZ :  2.47/16 11/09/2003  15.10.02  by  Michael Scheer
*CMZ :  2.34/09 18/09/2001  22.52.09  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.38.34  by  Michael Scheer
*CMZ :  2.17/00 03/11/2000  09.47.59  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  12.09.17  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.16/00 27/05/2000  14.03.56  by  Michael Scheer
*CMZ :  2.15/00 08/05/2000  13.32.10  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  2.13/02 14/12/99  16.24.13  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.47.16  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  12.37.49  by  Michael Scheer
*CMZ : 00.02/00 19/11/96  14.57.13  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  17.29.50  by  Michael Scheer
*CMZ : 00.01/06 01/02/95  16.35.43  by  Michael Scheer
*CMZ : 00.01/04 29/11/94  10.17.51  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.07.25  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.56  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.26  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PININ
*KEEP,GPLHINT.
*KEND.

*KEEP,SOURCEF90U.

      USE SOURCEF90

*KEEP,OBSERVF90U.

	USE OBSERVF90
*KEEP,WFOLDF90U.

	USE WFOLDF90
*KEND.

C--- INITIALIZE GRID OF OBERSERVATION POINTS OF PINHOLE

      IMPLICIT NONE

*KEEP,CMPARA.

      INTEGER LIDIMP,NOMDIMP,IBFDIM4P,IBFDIM2P,NXPANP,NYPANP
      INTEGER NBTABP,NDMASHZP,NDMASHYP,MAXFOUR,NFOURD

      INTEGER NWMAXP,NDFREQP,NDOBSVP,NDARGUP,NDWSOUP
      INTEGER NGCOEFP,NDOBSVZP,NDOBSVYP,NSPLINEP,NDPOLP,NDSPARP

      INTEGER N3DPOIP,N2DHPOIP,NPHPOIP


      PARAMETER (LIDIMP=501)    ! MAX. NUMBER OF SOURCES
      PARAMETER (NDPOLP=101)    ! MAX. NUMBER OF POLES
                                ! IN SR BEAMPOW POLES MAY BE SPLITTED, I.E.
                                ! NUMBER OF POLES ARE NOT NUMBER OF
                                ! HARDWAREPOLES
      !
      PARAMETER (NWMAXP=2**15)    !NUMBER OF POINTS OF REFERENCE ORBIT
      PARAMETER (NDFREQP=1100001)    !NUMBER OF FREQUENCES FOR SPECTRUM-CALCULATION
      PARAMETER (NDSPARP=1024)   !TO HAVE SOME SPARE IN BUFFER
      PARAMETER (NDWSOUP=NDSPARP+1024)!NUMBER OF STEPS PER SOURCE OF
                                      !(MINI-TRAJEKTORY)
               !SEE ALSO PARAMETER NBADDP IN SOURCE.CMN
               !FOR INTEGRATION BUFFER IN SR SPECTRUM
      !
      PARAMETER (NDARGUP=NDSPARP+1024)!NUMBER OF INTEGRATION STEPS PER SOURCE
                                      !MUST BE LOWER OR EQUAL NDWSOUP
      !
      PARAMETER (NDOBSVZP=2048)    !NUMBER OF OBSERVATION POINTS IN Z-DIRECTION
      ! SAME FOR Z AND Y, OR PROBLEMS...
      PARAMETER (NDOBSVYP=NDOBSVZP)    !NUMBER OF OBSERVATION POINTS IN Y-DIRECTION
      PARAMETER (NDOBSVP=NDOBSVZP*NDOBSVYP) !NUMBER OF OBSERVATION POINTS
      PARAMETER (NSPLINEP=NDOBSVZP+NDOBSVYP)!FOR SR FSPLINE

      PARAMETER (NGCOEFP=32)         !NUMBER FOURIER COEFFICIENT FOR GAUSSIAN FOR
      !                                !FOLDING PROCEDURE
      PARAMETER (NDMASHZP=NDOBSVZP/2)   !NUMBER OF ADJACENT MASHES FOR FOLDING
      PARAMETER (NDMASHYP=NDOBSVYP/2)   !NUMBER OF ADJACENT MASHES FOR FOLDING

      PARAMETER(MAXFOUR=2**15)  !NUMBER OF FOURIER COEFFICIENTS FOR B-FIELD
      PARAMETER (NFOURD=2*MAXFOUR) !NUMBER OF FOURIER COEFFICIENTS FOR SR FOURWLS

      INTEGER NBMESSXP,NBMESSYP,NBMESSZP
      PARAMETER (NBMESSXP=3,NBMESSYP=3,NBMESSZP=3)

      PARAMETER(N3DPOIP=NBMESSXP*NBMESSYP*NBMESSZP) !NUMBER DATA POINTS FOR BPOLY3D
      PARAMETER(N2DHPOIP=NBMESSXP*NBMESSYP*NBMESSZP) !NUMBER DATA POINTS FOR BPOLY3D
      PARAMETER(NPHPOIP=NBMESSXP*NBMESSYP*NBMESSZP) !NUMBER DATA POINTS FOR BPOLY3D

      INTEGER NLIOBFRP,NOBFRP,NLIFRP,NLIOBZP,NLIOBP
      INTEGER ILIOB,ILIOBFR,IOBFR,ILIOBZ,ILIFR,IFROB
        PARAMETER (NLIOBFRP=600000)  !GENERAL DIMENSION OF SPECTRUM ARRAYS
      PARAMETER (NLIOBP=NLIOBFRP)
      PARAMETER (NLIOBZP=NLIOBFRP)
      PARAMETER (NOBFRP=NLIOBFRP)
      PARAMETER (NLIFRP=NLIOBFRP)

C------------ OLD STUFF (??) --------------------------------

      PARAMETER (NOMDIMP=1)   ! MAX. ANZAHL VON FREQUENZEN IM SPEKTRUM
      PARAMETER (IBFDIM4P=1)    ! BUFFER GROESSE FUER RBUFF4,IBUFF4
      PARAMETER (IBFDIM2P=8)    ! BUFFER GROESSE FUER RBUFF2
      PARAMETER (NXPANP=1)      ! PANDIRA X-DIMENSION
      PARAMETER (NYPANP=1)      ! PANDIRA Y-DIMENSION
      PARAMETER (NBTABP=2**14+2**13)   ! BTAB DIMENSION

c20180509      REAL*4 RBUFF4(4,IBFDIM4P)
c      REAL*4 RBUFF2(2,IBFDIM2P)
c      INTEGER IBUFF4(4,IBFDIM4P)
c      INTEGER IBUFF2(2,IBFDIM2P)
c      EQUIVALENCE (RBUFF4(1,1),IBUFF4(1,1)),(RBUFF2(1,1),IBUFF2(1,1))
c      COMMON/BUFCOM/RBUFF4,RBUFF2,ILIOB,ILIOBFR,IOBFR,ILIOBZ,ILIFR,IFROB
*KEEP,CONTRL.
      DOUBLE PRECISION ECPHOTON,XSTART,XINTER,XSTOP,XIANF,XIEND,GWINFC,BL0CUT,BL0HYS
     &  ,DMYGAMMA,DMYBETA,VXIN,VYIN,DMYCUR,BANWID,D1MBETA
     &  ,VZIN,YSTART,ZSTART,DEVLEN,DEVLEN2,DMYENERGY,GMOM,EMOM
     &  ,XSTARTH,XSTOPH,DBRHO,BMOVECUT
     &  ,BXSTART,BYSTART,BZSTART,AXSTART,AYSTART,AZSTART
     &  ,BXSTOP,BYSTOP,BZSTOP,AXSTOP,AYSTOP,AZSTOP,XBSYM

      INTEGER       LIDIM,NBDIM,NOMDIM,IBFDIM4,IBFDIM2,IANZPL,LUNGFI,
     &  LUNGFO,IWFILSP0,IWFILSPF,IWFILPOW,
     &  NXPAN,NYPAN,
     &  ICODE,KOUT,KBFELD,KHALBA,KHALBASY,KUNDUGAP,KBEXTERN,KBAMWLS,
     &  IRFILF,IBGAUSS,KCIRC,
     &  KSIGN,ICHECK,KBETAX,IFORM0,IFORM,
     &  IWFILB0,IRFILB0,IWFILF,IWSECTMAGS,IRFILP,IRFILB,ITRAKT,
     &  IWFILT0,IRFILT0,IWFILL0,IRFILL0,
     &  IEXPL0,NLPOI,IWFILL,IWFILA,IRFILA,ISAVLO,
     &  IOPTIC,IGENFUN,IDISPER,IEMIT,IERZANA,IEMIAHW,
     &  IBHARD,IERZFUN,IRANDO,IBHTRACK,IBHELM,KBGENESIS,ISPLINE,
     &  IBSYM,IWLSOPT,IKBFORM,IRBTAB,NBTAB,IJUST,MYINUM,NSTEPMX,
     &  ICTEST ,IHFOLD,IBSYMY,IBSYMZ,
     &  NWMAX,ISPEC,IVELOFIELD,NDFREQ,NDOBSV,IFREQ2P,IWFILINT,JWFILINT,
     &  IHPIN,
     &  NDPAWC,IPIN,IPOLA,ifold,ibunch,ISPECANA,NGCOEF,ISPECMODE,
     &  NDMASHZ,NDMASHY,NDOBSVZ,NDOBSVY,ISIGUSR,ISPECINT,
     &  IPINALL,IHBOOK,IHINDEX,
     &  IHTRACK,IHTRSMP,IHTRACKM,IHBETA,IUNIT,IUNITS,IF1DIM,IHFREQ,
     &  IPOWER,NDPOL,IBEAMPOL,KUCROSS,KMAGSEQ,IMGSQF,KMAGCOR,
     &  IWBTAB,IABEND,irbmap6,iwbmap,
     &  IWFILFL0,IWFILFLF,IRPHI,
     &  IUSEM,IPINCIRC,ISTOKES,irbtabzy,ifourbtabzy,IRBTABXYZ,ISPECANAF,IBRILL,
     &  IWFILSTO,IWFLSTOF,IWFLSTOE,IWFLSTOEF,IRFILSP0,IRFILSTO,ISPECSUM,
     &  IWFILS,IWFILSF,IWFILSE,IWFILSEF,
     &  IWFILB,IWFILBF,IWFILBE,IWFILBEF,
     &  IWFILBRILL,IWFILBRILLF,IWFILBRILLE,IWFILBRILLEF,
     &  KELLIP,IWFILRAY,KELLANA,
     &  IPHASE,IDOSE,IPHOTON,IEFOLD,IDESYNC,ISPECDIP,
     &  IUOUT,IUNAME,IUSTEP,IPHASEANA,
     &  KBPOLYH,IGFLOAT
     &  ,KBPOLY3D,IWBPOLY3D
     &  ,KBPOLY2DH,IWBPOLY2DH
     &  ,KBPharm,IWBpharm
     &  ,KBREC,KBPOLYMAG,IBATCH,kbundumag,kbundumap,kbunduverb
     &  ,IMAGSPLN,IHINPUT,IHOUTP,IMHBCOM,IMAGJOB
     &  ,IW_BLEN,IW_BLENF,IW_CIRC
     &  ,IHLIMIT_C,IHISINI_C,ihisascii,iroottrees,IBFORCE
     &  ,iampli,kampli,IAMPJIT,IUNDULATOR,IWIGGLER,IEXPERT,IBSUPER,IBERROR
     &  ,ieneloss,iefield,iroothdf5,nocern,iadjust,iwarnbmap,iwarnmyb,mthreads

      CHARACTER(65) CODE
      CHARACTER     CHISASCII

      COMMON/CONTRL/ECPHOTON,
     &  XSTART,XINTER,XSTOP,XIANF,XIEND,GWINFC,BL0CUT,BL0HYS
     &  ,DMYGAMMA,DMYBETA,VXIN,VYIN,DMYCUR,BANWID,D1MBETA
     &  ,VZIN,YSTART,ZSTART,DEVLEN,DEVLEN2,DMYENERGY,GMOM,EMOM
     &  ,XSTARTH,XSTOPH,DBRHO,BMOVECUT
     &  ,BXSTART,BYSTART,BZSTART,AXSTART,AYSTART,AZSTART
     &  ,BXSTOP,BYSTOP,BZSTOP,AXSTOP,AYSTOP,AZSTOP,XBSYM
     &  ,LIDIM,NBDIM,NOMDIM,IBFDIM4,IBFDIM2,IANZPL,LUNGFI,
     &  LUNGFO,IWFILSP0,IWFILSPF,IWFILPOW,
     &  NXPAN,NYPAN,
     &  ICODE,KOUT,KBFELD,KHALBA,KHALBASY,KUNDUGAP,KBEXTERN,KBAMWLS,
     &  IRFILF,IBGAUSS,KCIRC,
     &  KSIGN,ICHECK,KBETAX,IFORM0,IFORM,
     &  IWFILB0,IRFILB0,IWFILF,IWSECTMAGS,IRFILP,IRFILB,ITRAKT,
     &  IWFILT0,IRFILT0,IWFILL0,IRFILL0,
     &  IEXPL0,NLPOI,IWFILL,IWFILA,IRFILA,ISAVLO,
     &  IOPTIC,IGENFUN,IDISPER,IEMIT,IERZANA,IEMIAHW,
     &  IBHARD,IERZFUN,IRANDO,IBHTRACK,IBHELM,KBGENESIS,ISPLINE,
     &  IBSYM,IWLSOPT,IKBFORM,IRBTAB,NBTAB,IJUST,MYINUM,NSTEPMX,
     &  ICTEST ,IHFOLD,IBSYMY,IBSYMZ,
     &  NWMAX,ISPEC,IVELOFIELD,NDFREQ,NDOBSV,IFREQ2P,IWFILINT,JWFILINT,
     &  IHPIN,
     &  NDPAWC,IPIN,IPOLA,ifold,ibunch,ISPECANA,NGCOEF,ISPECMODE,
     &  NDMASHZ,NDMASHY,NDOBSVZ,NDOBSVY,ISIGUSR,ISPECINT,
     &  IPINALL,IHBOOK,IHINDEX,
     &  IHTRACK,IHTRSMP,IHTRACKM,IHBETA,IUNIT,IUNITS,IF1DIM,IHFREQ,
     &  IPOWER,NDPOL,IBEAMPOL,KUCROSS,KMAGSEQ,KMAGCOR,IMGSQF,
     &  IWBTAB,IABEND,irbmap6,iwbmap,
     &  IWFILFL0,IWFILFLF,IRPHI,
     &  IUSEM,IPINCIRC,ISTOKES,irbtabzy,ifourbtabzy,IRBTABXYZ,ISPECANAF,IBRILL,
     &  IWFILSTO,IWFLSTOF,IWFLSTOE,IWFLSTOEF,IRFILSP0,IRFILSTO,ISPECSUM,
     &  KELLIP,IWFILRAY,KELLANA,
     &  IPHASE,IDOSE,IPHOTON,IEFOLD,IDESYNC,ISPECDIP,IPHASEANA,
     &  IUOUT,IUNAME,IUSTEP,
     &  KBPOLYH,IGFLOAT
     &  ,KBPOLY3D,IWBPOLY3D
     &  ,KBPOLY2DH,IWBPOLY2DH,
     &  KBPharm,IWBpharm,
     &  IWFILS,IWFILSF,IWFILSE,IWFILSEF,
     &  IWFILB,IWFILBF,IWFILBE,IWFILBEF,
     &  IWFILBRILL,IWFILBRILLF,IWFILBRILLE,IWFILBRILLEF
     &  ,KBREC,KBPOLYMAG,IBATCH,kbundumag,kbundumap,kbunduverb
     &  ,IMAGSPLN,IHINPUT,IHOUTP,IMHBCOM, IMAGJOB
     &  ,IW_BLEN,IW_BLENF,IW_CIRC
     &  ,IHLIMIT_C,IHISINI_C,ihisascii,iroottrees,IBFORCE
     &  ,iampli,kampli,IAMPJIT,IUNDULATOR,IWIGGLER,IEXPERT,IBSUPER,IBERROR
     &  ,ieneloss,iefield,iroothdf5,iwarnbmap,iwarnmyb,mthreads
     &  ,CODE,CHISASCII


      NAMELIST/CONTRL/CODE,ECPHOTON,MYINUM,NSTEPMX,DMYENERGY,VXIN,VYIN,VZIN,BMOVECUT
     &  ,KBFELD,KHALBA,KHALBASY,KUNDUGAP,IRFILF,IRBTAB,
     &  KBEXTERN,KBAMWLS,IBGAUSS,KCIRC,IRPHI,
     &  xstart,XINTER,xstop,xianf,xiend,
     &  KSIGN,ICHECK,KBETAX,IFORM0,IFORM,
     &  IWFILB0,IRFILB0,IWFILF,IWSECTMAGS,IRFILP,IRFILB,ITRAKT,
     &  IWFILT0,IRFILT0,IWFILL0,GWINFC,BL0CUT,BL0HYS,
     &  IRFILL0,
     &  IEXPL0,NLPOI,IWFILL,IWFILA,IRFILA,ISAVLO,IOPTIC,
     &  IBHARD,IERZFUN,IERZANA,IGENFUN,
     &  IRANDO,IBHTRACK,IBHELM,KBGENESIS,ISPLINE,IDISPER,
     &  IBSYM,XBSYM,IWLSOPT,IKBFORM,IJUST,IBSYMY,IBSYMZ
     &  ,YSTART,ZSTART,DMYCUR,BANWID
     &  ,ICTEST,IEMIT,IEMIAHW,IHFOLD
     &  ,ISPEC,IVELOFIELD,IWFILSP0,IFREQ2P,IWFILINT,JWFILINT,IWFILSPF
     &  ,IWFILPOW
     &  ,IHPIN,IPOLA,ifold,ibunch,IPIN,ISPECANA,ISIGUSR,ISPECINT,ISPECMODE
     &  ,IPINALL,IHINDEX
     &  ,IHTRACK,IHTRSMP,IHTRACKM,IHBETA,IUNIT,IUNITS,IF1DIM,IHFREQ
     &  ,XSTARTH,XSTOPH,IPOWER,IBEAMPOL,KUCROSS,KMAGSEQ,KMAGCOR,IMGSQF,
     &  IWBTAB,irbmap6,iwbmap
     &  ,IWFILFL0,IWFILFLF
     &  ,IUSEM,IPINCIRC,ISTOKES,irbtabzy,ifourbtabzy,IRBTABXYZ,ISPECANAF,IBRILL
     &  ,IWFILSTO,IWFLSTOF,IWFLSTOE,IWFLSTOEF,IRFILSP0,IRFILSTO,ISPECSUM
     &  ,KELLIP,IWFILRAY,KELLANA
     &  ,IPHASE,IDOSE,IPHOTON,IEFOLD,IDESYNC,ISPECDIP,IPHASEANA
     &  ,IUOUT,IUNAME,IUSTEP
     &  ,KBPOLYH
     &  ,KBPOLY3D,IWBPOLY3D
     &  ,KBPOLY2DH,IWBPOLY2DH,
     &  KBPharm,IWBpharm,
     &  IWFILS,IWFILSF,IWFILSE,IWFILSEF,
     &  IWFILB,IWFILBF,IWFILBE,IWFILBEF,
     &  IWFILBRILL,IWFILBRILLF,IWFILBRILLE,IWFILBRILLEF
     &  ,KBREC,KBPOLYMAG,kbundumag,kbundumap,kbunduverb
     &  ,IMAGSPLN,IHINPUT,IHOUTP,IMHBCOM, ihisascii,iroottrees,
     &  ieneloss,iefield,iroothdf5,
     &  iampli,kampli,IAMPJIT,IUNDULATOR,IWIGGLER,IEXPERT,IBSUPER,IBERROR,ibforce,
     &  nocern,iadjust,mthreads,
     &  CHISASCII
*KEEP,MYFILES.
      CHARACTER(128) FILEB0,FILET0,FILEL0,FILEL,FILEA,FILEP,FILEO,FILECOD,
     &  FILEH,FILED,FILEF,FILETB,FILEJ,FILETR,FILEI,FILEBE
     &  ,FILEFR,FILEOB,FILEGFO,FILESP0,FILEINT,FILEHB
     &  ,FILEWB,FILEZZPYYP,FILESPF,FILEWBT,FILEPOW
     &  ,FILEABS,FILEMG,FILEFF
     &  ,FILEFL0,FILEFLF,filetbx
     &  ,FILETBZ,FILESTO,FILESTOF,FILESTOE,FILESTOEF,FILERAY
     &  ,FILEAM,FILEAMO
     &  ,FILES,FILESE,FILESF,FILESEF
     &  ,FILEC,FILECE,FILECF,FILECEF
     &  ,FILEBRILL,FILEBRILLE,FILEBRILLF,FILEBRILLEF
     &  ,FILEREC,FILEBMAP
     &  ,FILE3DCOE,FILE3DFIT
     &  ,FILE2DHFIT
     &  ,FILEPHFIT
     &  ,FILEPH,FILEAMPLI,FILEFTH,FILEFTV,FILEGENL,FILEGENI

C5.10.95      REAL*4 B(NBDIMP,3,3,3),RB(NBDIMP,3)
      REAL*4 DYB0,DZB0

      INTEGER LUNB0,NXB0,NYB0,NZB0,LUNT0,LUNL0,LUNL,LUNA,LUNP,LUNO,
     &  LUND,LUNCOD,LUNH,LUNF,LUNTB,LUNJ,LUNTR,LUNBE
     &  ,LUNOB,LUNFR,LUNSP0,LUNINT,LUNHB,LUNWB,LUNZZPYYP
     &  ,LUNSPF,LUNWBT
     &  ,LUNABS,LUNMG,LUNEFF
     &  ,LUNFL0,LUNFLF
     &  ,LUNTBZ,LUNSTO,LUNSTOF,LUNRAY
     &  ,LUNAM,LUNAMO
     &  ,LUNS,LUNSE,LUNSF,LUNSEF
     &  ,LUNC,LUNCE,LUNCF,LUNCEF,LUNREC,LUNBMAP
     &  ,LUN3DCOE,LUN3DFIT
     &  ,LUN2DHFIT
     &  ,LUNPHFIT
     &  ,LUNPH,LUNAMPLI,LUNFT,LUNPOW,LUNGENL,LUNGENI

      COMMON/MYFILES/FILEB0,LUNB0,DYB0,DZB0,NXB0,NYB0,NZB0,
     &  FILET0,LUNT0,FILEL0,LUNL0,
     &  FILEL,LUNL,FILEA,LUNA,
     &  FILEP,LUNP,FILEO,LUNO,FILED,LUND,
     &  FILECOD,LUNCOD,
     &  FILEH,LUNH,
     &  FILEF,LUNF,
     &  FILETB,LUNTB
     &  ,FILEJ,LUNJ
     &  ,FILETR,LUNTR,FILEI,FILEGFO,FILEBE,LUNBE
     &  ,FILEOB,LUNOB
     &  ,FILEFR,LUNFR
     &  ,FILESP0,LUNSP0
     &  ,FILEPOW,LUNPOW
     &  ,FILESPF,LUNSPF
     &  ,FILEINT,LUNINT
     &  ,FILEHB,LUNHB
     &  ,FILEWB,LUNWB
     &  ,FILEZZPYYP,LUNZZPYYP
     &  ,FILEWBT,LUNWBT
     &  ,FILEABS,LUNABS,FILEFF,LUNEFF
     &  ,FILEMG,LUNMG
     &  ,FILEFL0,FILEFLF
     &  ,LUNFL0,LUNFLF
     &  ,FILETBZ,LUNTBZ,filetbx
     &  ,LUNSTO,FILESTO
     &  ,LUNSTOF,FILESTOF,FILESTOE,FILESTOEF
     &  ,LUNRAY,FILERAY
     &  ,LUNAM,FILEAM
     &  ,LUNAMO,FILEAMO
     &  ,LUNS,FILES
     &  ,LUNSE,FILESE
     &  ,LUNSF,FILESF
     &  ,LUNSEF,FILESEF
     &  ,LUNC,FILEC
     &  ,LUNCE,FILECE
     &  ,LUNCF,FILECF
     &  ,LUNCEF,FILECEF
     &  ,FILEBRILL,FILEBRILLE,FILEBRILLF,FILEBRILLEF
     &  ,LUNREC,FILEREC
     &  ,LUNBMAP,FILEBMAP
     &  ,LUN3DCOE,LUN3DFIT
     &  ,FILE3DCOE,FILE3DFIT
     &  ,LUN2DHFIT,FILE2DHFIT
     &  ,LUNPHFIT,FILEPHFIT
     &  ,LUNPH,FILEPH,LUNGENL,LUNGENI,FILEGENL,FILEGENI
     &  ,LUNAMPLI,FILEAMPLI,LUNFT,FILEFTH,FILEFTV
C5.10.95     &                 B,RB,

      NAMELIST/MYFILES/FILEB0,LUNB0,DYB0,DZB0,NYB0,NZB0,NXB0,
     &  FILET0,LUNT0,FILEL0,LUNL0,
     &  FILEL,LUNL,FILEA,LUNA,FILEP,LUNP,
     &  FILEO,LUNO,FILED,LUND,
     &  FILECOD,LUNCOD,
     &  FILEH,LUNH,
     &  FILEF,LUNF,
     &  FILETB,LUNTB
     &  ,FILEJ,LUNJ
     &  ,FILETR,LUNTR,FILEBE,LUNBE
     &  ,FILEOB,LUNOB
     &  ,FILEFR,LUNFR
     &  ,FILESP0,LUNSP0
     &  ,FILEPOW,LUNPOW
     &  ,FILEINT,LUNINT
     &  ,FILEHB,LUNHB
     &  ,FILEWB,LUNWB
     &  ,FILEZZPYYP,LUNZZPYYP
     &  ,FILESPF,LUNSPF
     &  ,FILEWBT,LUNWBT
     &  ,FILEABS,LUNABS,FILEFF,LUNEFF
     &  ,FILEMG,LUNMG
     &  ,FILEFL0,FILEFLF
     &  ,LUNFL0,LUNFLF
     &  ,FILETBZ,LUNTBZ,filetbx
     &  ,LUNSTO,FILESTO
     &  ,LUNSTOF,FILESTOF,FILESTOE,FILESTOEF
     &  ,LUNRAY,FILERAY
     &  ,LUNAM,FILEAM
     &  ,LUNAMO,FILEAMO
     &  ,LUNS,FILES
     &  ,LUNSE,FILESE
     &  ,LUNSF,FILESF
     &  ,LUNSEF,FILESEF
     &  ,LUNC,FILEC
     &  ,LUNCE,FILECE
     &  ,LUNCF,FILECF
     &  ,LUNCEF,FILECEF
     &  ,FILEBRILL,FILEBRILLE,FILEBRILLF,FILEBRILLEF
     &  ,LUNREC,FILEREC
     &  ,LUNBMAP,FILEBMAP
     &  ,LUN3DCOE,LUN3DFIT
     &  ,FILE3DCOE,FILE3DFIT
     &  ,LUN2DHFIT,FILE2DHFIT
     &  ,LUNPHFIT,FILEPHFIT
     &  ,LUNPH,FILEPH,LUNGENL,LUNGENI,FILEGENL,FILEGENI
     &  ,LUNAMPLI,FILEAMPLI,LUNFT,FILEFTH,FILEFTV
*KEEP,OBSERVF90.

      DOUBLE PRECISION PINCEN(3),PINW,PINH,pinr,pinrad,
     &  OBSVDZ,OBSVDY,OBS1X,OBS1Y,OBS1Z,OBSVDR,OBSVDPHI
     &  ,pinhsc,pinwsc,rpinsph

      INTEGER NOBSV,NOBSVZ,NOBSVY,MPINZ,MPINY,MPINR,MPINPHI
     &  ,mpinyorig,mpinzorig
      INTEGER MOBSV,MOBSVZ,MOBSVY,MEDGEZ,MEDGEY
      INTEGER MMEDGEZ,MMEDGEY,IRFILOB,MEDGER
     &  ,IPBRILL,ICBRILL
     &  ,NOBSVRPHI,NOBSVR,NOBSVPHI
     &  ,MOBSVRPHI,MOBSVR,MOBSVPHI,IQUADPHI

      COMMON/OBSERV/OBSVDZ,OBSVDY,OBS1X,OBS1Y,OBS1Z
     &  ,PINCEN,PINW,PINH,pinr,pinrad,OBSVDR,OBSVDPHI
     &  ,pinhsc,pinwsc
     &  ,rpinsph
     &  ,NOBSV,NOBSVZ,NOBSVY,MPINZ,MPINY,MPINR,MPINPHI
     &  ,mpinyorig,mpinzorig
     &  ,MOBSV,MOBSVZ,MOBSVY,MEDGEZ,MEDGEY
     &  ,MMEDGEZ,MMEDGEY,IRFILOB,MEDGER
     &  ,NOBSVRPHI,NOBSVR,NOBSVPHI
     &  ,MOBSVRPHI,MOBSVR,MOBSVPHI
     &  ,IPBRILL,ICBRILL,IQUADPHI

      NAMELIST/PINHOLE/PINCEN,MPINZ,MPINY,MPINR,MPINPHI,OBSVDR,OBSVDPHI,
     &  PINW,PINH,pinr,pinrad,MEDGEZ,MEDGEY
     &  ,MMEDGEZ,MMEDGEY,IRFILOB,OBS1X,OBS1Y,OBS1Z,OBSVDZ,OBSVDY
     &  ,IPBRILL,IQUADPHI
     &  ,pinhsc,pinwsc
     &  ,rpinsph
*KEEP,DEPOLA.
      DOUBLE PRECISION TAUPOL01G
      DOUBLE PRECISION RDIPOL,BDIPOL,UMFANG,TAUPOL,POLFAC,TAUKRIT
      DOUBLE PRECISION POLLVEP,TAUPLEP
      DOUBLE PRECISION BETFUN,BETFUNV,DI2RING,DI5RING,BAHWMX,BAHWMN,DBAHW,
     &  EMIKRIT,DISP0,DDISP0,DI4RING,X0BET,POLLEV,xbetfun,EPS0H,EPS0V,phaseh,
     &  phasev,BETAH,BETAPH,BETAV,BETAPV,DELGAM,DI3RING,DI3WLS,SIGEE

      integer ibetback

      COMMON/DEPOLAC/RDIPOL,BDIPOL,UMFANG,TAUPOL,POLFAC,TAUKRIT,
     &  BETFUN,BETFUNV,DI2RING,DI5RING,BAHWMX,BAHWMN,DBAHW,
     &  TAUPLEP,POLLVEP,EMIKRIT,DISP0,DDISP0,DI4RING,X0BET,POLLEV,EPS0H,EPS0V,
     &  phaseh,phasev,BETAH,BETAPH,BETAV,BETAPV,DELGAM,xbetfun,TAUPOL01G,
     &  DI3RING,DI3WLS,SIGEE,ibetback

      NAMELIST/DEPOLA/RDIPOL,UMFANG,TAUKRIT,BETFUN,BETFUNV,X0BET,
     &  DI2RING,DI5RING ,DISP0,DDISP0,DI4RING,EPS0H,EPS0V,phaseh,phasev,
     &  xbetfun,BETAH,BETAPH,BETAV,BETAPV,DELGAM,ibetback
*KEEP,WFOLDF90.
C      INTEGER max

      INTEGER NGFOURZ,NGFOURY,ISIGSTO
     &         ,NSIGE

      DOUBLE PRECISION GCOEFH(NGCOEFP,LIDIMP),GCOEFV(NGCOEFP,LIDIMP)
      DOUBLE PRECISION XKGAUSS(NGCOEFP,LIDIMP),YKGAUSS(NGCOEFP,LIDIMP)

      DOUBLE PRECISION WSIG2Z(max(10,LIDIMP)),WSIG2Y(max(10,LIDIMP))
      DOUBLE PRECISION DSIG2Z(max(10,LIDIMP)),DSIG2Y(max(10,LIDIMP))
     &        ,WSIGZ(max(10,LIDIMP)),WSIGY(max(10,LIDIMP))
     &        ,BSIGZ(max(10,LIDIMP)),BSIGY(max(10,LIDIMP))
     &        ,SIGRC,SIGRPC
     &        ,BSIGZP(max(10,LIDIMP)),BSIGYP(max(10,LIDIMP))
     &        ,DSIGZ(max(10,LIDIMP)),DSIGY(max(10,LIDIMP))
     &        ,DGSIGZ(max(10,LIDIMP)),DGSIGY(max(10,LIDIMP))

      DOUBLE PRECISION ESPREAD


      COMMON/WFOLDC/
     &         GCOEFH,GCOEFV
     &        ,XKGAUSS,YKGAUSS
     &        ,WSIG2Z,WSIG2Y
     &        ,DSIG2Z,DSIG2Y
     &        ,WSIGZ,WSIGY
     &        ,BSIGZ,BSIGY
     &        ,SIGRC,SIGRPC
     &        ,BSIGZP,BSIGYP
     &        ,DSIGZ,DSIGY
     &        ,DGSIGZ,DGSIGY
     &        ,ESPREAD
     &        ,NGFOURZ,NGFOURY,ISIGSTO
     &        ,NSIGE


      NAMELIST/WFOLDN/BSIGZ,BSIGY,BSIGZP,BSIGYP
     &        ,WSIGZ,WSIGY,DGSIGZ,DGSIGY,NGFOURZ,NGFOURY,ISIGSTO
     &        ,SIGRC,SIGRPC
     &        ,ESPREAD,NSIGE

*KEEP,SOURCEF90.

      INTEGER NBADDP,NROIP
      PARAMETER (NBADDP=100,NROIP=512)
      INTEGER NROI,NROIA,NIDSOURCE,NSOURCE,NDWSOU,NLPOIO,
     &          ISOURO,NSADD,NBUFF,NDSPAR,IROIASYEXP(NROIP)

      DOUBLE PRECISION ROIX(NROIP),ROIP(NROIP)

      COMMON/SOURCEC/ROIX,ROIP
     &,NROI,NROIA,NIDSOURCE,NSOURCE,NDWSOU,NLPOIO,ISOURO,NSADD,NBUFF,NDSPAR
     &,IROIASYEXP

      NAMELIST /ROIN/ NROI,ROIX,ROIP,IROIASYEXP
*KEND.


      INTEGER IOB,IY,IZ,INCZ,INCY,ISOUR,INCZMX,INCYMX
      INTEGER N2POWZ,N2POWY,ICDUMZ,ICDUMY

      DOUBLE PRECISION OBSVDUM(3,1),x,y,z,xn,yn,zn,r,r0,pinwo,pinho,pinro

      INTEGER IFPHIR_A,iobsv
      DATA IFPHIR_A/0/

      pinwo=pinw
      pinho=pinh
      pinro=pinr

      IF (IPINCIRC.NE.0) THEN
        PINW=2.D0*PINR
        PINH=2.D0*PINR
      ENDIF !IPINCIRC

      IF (OBSVDZ.NE.0.D0) THEN
        MPINZ=NINT(PINW/OBSVDZ)+1
        if (mpinz.lt.3) then
          write(6,*)"*** Warning in pinin: OBSVDZ > PINW/2., will be adjusted ***"
          write(lungfo,*)"*** Warning in pinin: OBSVDZ > PINW/2., will be adjusted ***"
          mpinz=3
          obsvdz=pinw/2.0d0
        endif
        PINW=OBSVDZ*(MPINZ-1)
      ENDIF !(OBSVDZ.NE.0.D0)

      IF (OBSVDY.NE.0.D0) THEN
        MPINY=NINT(PINH/OBSVDY)+1
        if (mpiny.lt.3) then
          write(6,*)"*** Warning in pinin: OBSVDY > PINH/2., will be adjusted ***"
          write(lungfo,*)"*** Warning in pinin: OBSVDY > PINH/2., will be adjusted ***"
          mpiny=3
          obsvdy=pinh/2.0d0
        endif
        PINH=OBSVDY*(MPINY-1)
      ENDIF !(OBSVDZ.NE.0.D0)

      IF (IPINCIRC.NE.0) THEN
        PINR=min(pinw,pinh)/2.0d0
      ENDIF !IPINCIRC

      if (abs(pinwo-pinw).gt.pinw/10000.0d0) then
        print*,"*** Warning in PINI: PINW is adjusted, according to MPINZ and DOBSVDZ ***"
      endif

      if (abs(pinho-pinh).gt.pinh/10000.0d0) then
        print*,"*** Warning in PINI: PINH is adjusted, according to MPINY and DOBSVDY ***"
      endif

      if (abs(pinro-pinr).gt.pinr/10000.0d0) then
        print*,"*** Warning in PINI: PINR is adjusted, according to PINW and PINH ***"
      endif

      IF (PINCEN(2).EQ.9999.) THEN
        IF     (IPBRILL.EQ.0) THEN
          PINCEN(2)=0.0
        ELSE IF (IPBRILL.EQ.1) THEN
          PINCEN(2)=PINH/2.D0
        ELSE IF (IPBRILL.EQ.2) THEN
          PINCEN(2)=PINH/2.D0
        ELSE IF (IPBRILL.EQ.3) THEN
          PINCEN(2)=-PINH/2.D0
        ELSE IF (IPBRILL.EQ.4) THEN
          PINCEN(2)=-PINH/2.D0
        ENDIF !IPBRILL
      ELSE IF (PINCEN(2).EQ.-9999.) THEN
        PINCEN(2)=YSTART+VYIN/VXIN*(PINCEN(1)-XSTART)
      ENDIF !PINCEN(2)

      IF (PINCEN(3).EQ.9999.) THEN
        IF     (IPBRILL.EQ.0) THEN
          PINCEN(3)=0.0
        ELSE IF (IPBRILL.EQ.1) THEN
          PINCEN(3)=PINW/2.D0
        ELSE IF (IPBRILL.EQ.2) THEN
          PINCEN(3)=-PINW/2.D0
        ELSE IF (IPBRILL.EQ.3) THEN
          PINCEN(3)=-PINW/2.D0
        ELSE IF (IPBRILL.EQ.4) THEN
          PINCEN(3)=PINW/2.D0
        ENDIF !IPBRILL
      ELSE IF (PINCEN(3).EQ.-9999.) THEN
        PINCEN(3)=ZSTART+VZIN/VXIN*(PINCEN(1)-XSTART)
      ENDIF !PINCEN(3)

      IF(NDOBSVZ*NDOBSVY.NE.NDOBSV) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*) '*** ERROR IN PININ ***'
        WRITE(LUNGFO,*) 'DIMENSION DECLARATIONS NOT CONSISTENT'
        WRITE(LUNGFO,*) 'NDOBSV MUST BE EQUAL TO NDOBSVZ*NDOBSVY'
        WRITE(LUNGFO,*) 'CHANGE PARAMETER IN CMPARA.CMN'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
        WRITE(6,*) '*** ERROR IN PININ ***'
        WRITE(6,*) 'DIMENSION DECLARATIONS NOT CONSISTENT'
        WRITE(6,*) 'NDOBSV MUST BE EQUAL TO NDOBSVZ*NDOBSVY'
        WRITE(6,*) 'CHANGE PARAMETER IN CMPARA.CMN'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF

C--- DATA OF PINHOLE ARE TAKEN FORM NAMELIST

      IF (IF1DIM.NE.0.AND.
     &    (MEDGEZ.NE.0.OR.MMEDGEZ.NE.0.OR.MPINZ.NE.1)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININ ***'
        WRITE(LUNGFO,*)'FLAG IF1DIM SET BUT'
        WRITE(LUNGFO,*)
     &    'MEDGEZ OR MMEDGEZ IN NAMELIST PINHOLE NOT ZERO OR MPINZ NOT EQUAL ONE'
        WRITE(LUNGFO,*)'ADJUSTED TO APPROPRIATE VALUES'
        WRITE(LUNGFO,*)
c        WRITE(6,*)
c        WRITE(6,*)
c        WRITE(6,*)'*** WARNING IN PININ ***'
c        WRITE(6,*)'FLAG IF1DIM SET BUT'
c        WRITE(6,*)
c     &    'MEDGEZ OR MMEDGEZ IN NAMELIST PINHOLE NOT ZERO OR MPINZ NOT EQUAL ONE'
c        WRITE(6,*)'ADJUSTED TO APPROPRIATE VALUES'
c        WRITE(6,*)
        MPINZ=1
        MEDGEZ=0
        MMEDGEZ=0
      ENDIF

      IF (IUSEM.NE.0) THEN
C MAKE MPINZ AND MPINY EVEN
        MPINZ=(MPINZ+1)/2*2
        MPINY=(MPINY+1)/2*2
      ENDIF !IUSEM

      NOBSVY=MPINY
      NOBSVZ=MPINZ
      NOBSV=NOBSVY*NOBSVZ
      IF (NOBSV.LE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'*** NEGATIVE NUMBER OF OBSERVATION POINTS!'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'*** NEGATIVE NUMBER OF OBSERVATION POINTS!'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      OBSVDUM(1,1)=PINCEN(1)

      IF (MPINZ.GT.1) THEN
        OBSVDZ=PINW/DFLOAT(MPINZ-1)
        OBSVDUM(3,1)=PINCEN(3)-PINW/2.
      ELSE
        OBSVDZ=PINW
        OBSVDUM(3,1)=PINCEN(3)
      ENDIF

      IF (MPINY.GT.1) THEN
        OBSVDY=PINH/DFLOAT(MPINY-1)
        OBSVDUM(2,1)=PINCEN(2)-PINH/2.
      ELSE
        OBSVDY=PINH
        OBSVDUM(2,1)=PINCEN(2)
      ENDIF

      IF (IF1DIM.EQ.0.AND.(MEDGEZ.LT.1.OR.MEDGEY.LT.1)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'MEDGEZ OR MEDGEY IN NAMELIST PINHOLE LOWER THAN 1'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'MEDGEZ OR MEDGEY IN NAMELIST PINHOLE LOWER THAN 1'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF

      IF (MMEDGEZ.LT.0.OR.MMEDGEY.LT.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'MMEDGEZ OR MMEDGEY IN NAMELIST PINHOLE LOWER THAN 0'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'MMEDGEZ OR MMEDGEY IN NAMELIST PINHOLE LOWER THAN 0'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF

C- DETERMINE RMS VALUES FOR FOLDING

      IF (IFOLD.NE.0) CALL WSIGFOL

C- INCREASE PINHOLE

      MOBSVZ=NOBSVZ  !STORE VALUES
      MOBSVY=NOBSVY
      MOBSV=NOBSV

      INCZ=0
      INCY=0

      IF (IFOLD.NE.0.AND.IFOLD.NE.2) THEN

        INCZMX=-1
        INCYMX=-1
        DO ISOUR=1,NSOURCE
          INCZ=NINT(WSIGZ(ISOUR)*DGSIGZ(ISOUR)/OBSVDZ)
          IF(INCZ.GT.INCZMX) THEN
            INCZMX=INCZ
          ENDIF
          INCY=NINT(WSIGY(ISOUR)*DGSIGY(ISOUR)/OBSVDY)
          IF(INCY.GT.INCYMX) THEN
            INCYMX=INCY
          ENDIF
        ENDDO   !ISOUR

        INCZ=INCZMX
        INCY=INCYMX
        IF (IF1DIM.NE.0) INCZ=0

      ENDIF !IFOLD

C           FACTOR 2 FOR BOTH SIDES OF PINHOLE
      NOBSVZ=MOBSVZ+2*(INCZ+MEDGEZ+MMEDGEZ)
      NOBSVY=MOBSVY+2*(INCY+MEDGEY+MMEDGEY)

C--- ADJUST NUMBER OF POINTS ACCORDING TO POWERS OF IF FLAG IUSEM IS SET

      IF (IUSEM.NE.0) THEN

        N2POWZ=NINT(ALOG(FLOAT(NOBSVZ-1))/ALOG(2.))
        IF(NOBSVZ .GT. 2**N2POWZ) N2POWZ=N2POWZ+1
C150793        NOBSVZ=2**N2POWZ+1
        NOBSVZ=2**N2POWZ

        N2POWY=NINT(ALOG(FLOAT(NOBSVY-1))/ALOG(2.))
        IF(NOBSVY .GT. 2**N2POWY) N2POWY=N2POWY+1
C150793        NOBSVY=2**N2POWY+1
        NOBSVY=2**N2POWY

      ENDIF !IUSEM

      OBSVDUM(2,1)=OBSVDUM(2,1)-(NOBSVY-MOBSVY)/2*OBSVDY
      OBSVDUM(3,1)=OBSVDUM(3,1)-(NOBSVZ-MOBSVZ)/2*OBSVDZ

      NOBSV=NOBSVZ*NOBSVY
      IF (NOBSV.GT.4000000) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'*** MORE THAN 4 000 0000 OBSERVATION POINTS ***'
        WRITE(LUNGFO,*)'NUMBER OF POINTS IN Y:',NOBSVY
        WRITE(LUNGFO,*)'NUMBER OF POINTS IN Z:',NOBSVZ
        IF (IFOLD.NE.0) WRITE(LUNGFO,*)'PLEASE CHECK PINHOLE SIZE AND RELATION TO BEAM EMITTANCE'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'*** MORE THAN 1 000 000 OBSERVATION POINTS ***'
        WRITE(6,*)'NUMBER OF POINTS IN Y:',NOBSVY
        WRITE(6,*)'NUMBER OF POINTS IN Z:',NOBSVZ
        IF (IFOLD.NE.0) WRITE(6,*)'PLEASE CHECK PINHOLE SIZE AND RELATION TO BEAM EMITTANCE'
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      IF (IOBSV_A.NE.NOBSV) THEN
        IF (IOBSV_A.NE.0) DEALLOCATE(OBSV)
        ALLOCATE(OBSV(3,NOBSV))
        IOBSV_A=NOBSV
      ENDIF !(IOBSV_A.LT.NOBSV)

      IF (IOBSVZ_A.NE.NOBSVZ) THEN
        IF (IOBSVZ_A.NE.0) DEALLOCATE(OBSVZ)
        ALLOCATE(OBSVZ(NOBSVZ))
        IOBSVZ_A=NOBSVZ
      ENDIF !(IOBSVY_A.LT.NOBSVY)

      IF (IOBSVY_A.NE.NOBSVY) THEN
        IF (IOBSVY_A.NE.0) DEALLOCATE(OBSVY)
        ALLOCATE(OBSVY(NOBSVY))
        IOBSVY_A=NOBSVY
      ENDIF !(IOBSV_A.LT.NOBSV)

      OBSV(1,1)=OBSVDUM(1,1)
      OBSV(2,1)=OBSVDUM(2,1)
      OBSV(3,1)=OBSVDUM(3,1)

      IOB=0
      DO IY=1,NOBSVY
        DO IZ=1,NOBSVZ
          IOB=IOB+1
          OBSV(3,IOB)=OBSV(3,1)+OBSVDZ*(IZ-1)
          OBSV(2,IOB)=OBSV(2,1)+OBSVDY*(IY-1)
          OBSV(1,IOB)=OBSV(1,1)
        ENDDO
      ENDDO

      DO IZ=1,NOBSVZ
        OBSVZ(IZ)=OBSV(3,IZ)
      ENDDO

      DO IY=1,NOBSVY
        IOB=(IY-1)*NOBSVZ+1
        OBSVY(IY)=OBSV(2,IOB)
      ENDDO

      IF (IPIN.EQ.0) THEN
        ICBRILL=1
      ELSE !PIN
        IF (IPBRILL.EQ.0) THEN
          ICDUMZ=NOBSVZ/2+1
          ICDUMY=NOBSVY/2+1
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.1) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+1
          ICDUMY=(NOBSVY-MOBSVY)/2+1
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.2) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+MOBSVZ
          ICDUMY=(NOBSVY-MOBSVY)/2+1
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.3) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+1
          ICDUMY=(NOBSVY-MOBSVY)/2+MOBSVY
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.4) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+MOBSVZ
          ICDUMY=(NOBSVY-MOBSVY)/2+MOBSVY
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN PININ: IPBRILL WRONG ***'
          WRITE(LUNGFO,*)'*** CHECK NAMELIST PINHOLE ***'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN PININ: IPBRILL WRONG ***'
          WRITE(6,*)'*** CHECK NAMELIST PINHOLE ***'
          WRITE(6,*)
          WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
          STOP
        ENDIF
      ENDIF   !IPIN

      IF (MPINR.NE.0) IRPHI=1
      IF (IPINCIRC.NE.0.AND.IFPHIR_A.NE.NOBSV) THEN
        ALLOCATE(FPHIR(NOBSV))
        IF (IRPHI.NE.0) THEN
          ALLOCATE(XC(NOBSV))
          ALLOCATE(YC(NOBSV))
        ENDIF
        IFPHIR_A=NOBSV
      ENDIF !IPINCIRC

      IF (MPINR.NE.0.AND.IPIN.EQ.0) THEN
        MPINR=0
        WRITE(6,*)
        WRITE(6,*)
     &    '*** WARNING IN PININ: IPIN.EQ.0, THEREFORE MPINR SET ZERO ***'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING IN PININ: IPIN.EQ.0, THEREFORE MPINR SET ZERO ***'
        WRITE(LUNGFO,*)
      ENDIF

      IF (MPINR.NE.0.AND.IF1DIM.NE.0) THEN
        MPINR=0
        WRITE(6,*)
        WRITE(6,*)
     &    '*** WARNING IN PININ: IF1DIM.NE.0, THEREFORE MPINR SET ZERO ***'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING IN PININ: IF1DIM.EQ.0, THEREFORE MPINR SET ZERO ***'
        WRITE(LUNGFO,*)
      ENDIF

      IF (MPINR.NE.0) THEN
        CALL PININR
      ENDIF

      if (rpinsph.eq.-9999.0d0) rpinsph=pincen(1)

      if (rpinsph.ne.0) then
        r0=obsv(1,1)-rpinsph
        do iobsv=1,nobsv
          x=obsv(1,iobsv)-r0
          y=obsv(2,iobsv)
          z=obsv(3,iobsv)
          r=sqrt(x*x+y*y+z*z)
          xn=x/r
          yn=y/r
          zn=z/r
          obsv(1,iobsv)=r0+xn*rpinsph
          obsv(2,iobsv)=yn*rpinsph
          obsv(3,iobsv)=zn*rpinsph
        enddo
        if (mpinr.ne.0) then
          do iobsv=1,nobsvrphi
            x=obsvrphi(1,iobsv)-r0
            y=obsvrphi(2,iobsv)
            z=obsvrphi(3,iobsv)
            r=sqrt(x*x+y*y+z*z)
            xn=x/r
            yn=y/r
            zn=z/r
            obsvrphi(1,iobsv)=xn*rpinsph
            obsvrphi(2,iobsv)=yn*rpinsph
            obsvrphi(3,iobsv)=zn*rpinsph
          enddo
        endif

        DO iob=1,nobsv
          if (abs(obsv(2,iob)).lt.1.0d-12) obsv(2,iob)=0.0d0
          if (abs(obsv(3,iob)).lt.1.0d-12) obsv(3,iob)=0.0d0
        ENDDO

        DO IZ=1,NOBSVZ
          OBSVZ(IZ)=OBSV(3,IZ)
        ENDDO

        DO IY=1,NOBSVY
          IOB=(IY-1)*NOBSVZ+1
          OBSVY(IY)=OBSV(2,IOB)
        ENDDO

      endif

      RETURN
      END
