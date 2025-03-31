*CMZ :          05/09/2024  15.54.32  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  11.54.00  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      integer :: mthreads,ktime=1

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase_prop',1)

      if (ifieldprop_u.gt.0) then

        if (modepin_u.ne.1) then
          call urad_phase_prop_classic(mthreads)
        else
          call urad_phase_prop_mc(mthreads)
        endif

      else !(ifieldprop_u.gt.0)

        call urad_phase_prop_geo

      endif !(ifieldprop_u.gt.0)

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase_prop',0)

      return
      end
