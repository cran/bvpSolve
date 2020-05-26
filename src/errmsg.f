c ===================================================================================
c print R-messages
c ===================================================================================

      subroutine intpr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      integer data(ndata)

       if (ndata .eq. 1) then
         call rprintfi1(label//char(0), data(1))
       else if (ndata .eq. 2) then
         call rprintfi2(label//char(0), data(1), data(2))
       else if (ndata .eq. 3) then
         call rprintfi3(label//char(0), data(1), data(2), data(3))
       else if (ndata .ge. 4) then
         call rprintfi4(label//char(0),data(1),data(2),data(3),data(4))
       endif 
      end subroutine


      subroutine dblepr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      double precision data(ndata)

       if (ndata .eq. 1) then
         call rprintfd1(label//char(0), data(1))
       else if (ndata .eq. 2) then
         call rprintfd2(label//char(0), data(1), data(2))
       else if (ndata .eq. 3) then
         call rprintfd3(label//char(0), data(1), data(2), data(3))
       else if (ndata .ge. 4) then
         call rprintfd4(label//char(0),data(1),data(2),data(3),data(4))
       endif 
      end subroutine

      subroutine logpr(msg,len,ldum)
      integer len
      character (len=*) msg
      logical ldum
      character (len=8) logi

      if (ldum) then
         logi ='  TRUE '
      else
         logi ='  FALSE'
      endif

      call rprint(MSG//logi//char(0))

      end subroutine

c ===================================================================================
c print R-messages
c ===================================================================================
C just a string
      subroutine rprint(msg)
      character (len=*) msg

      call rprintf(msg//char(0))

      end subroutine 


C printing with one logical
      subroutine rprintl1(msg, l1)
      character (len=*) msg
      logical l1
      character (len=8) logi

       if (l1) then
         logi ='  TRUE '
       else
         logi ='  FALSE'
       endif
       call rprint(MSG//logi//char(0))

      end subroutine

C  two logicals
      subroutine rprintl2(msg, l1, l2)
      character (len=*) msg
      character (len=8) logi, logi2
      logical l1, l2

       if (l1) then
         logi ='  TRUE '
       else
         logi ='  FALSE'
       endif
      
       if (l2) then
         logi2 ='  TRUE '
       else
         logi2 ='  FALSE'
       endif
       call rprint(msg//logi//logi2//char(0))

      end subroutine 

C printing with one integer and a double
      subroutine rprintid(msg, i1, d1)
      character (len=*) msg
      double precision d1
      integer i1
         call rprintfd1(msg//char(0), d1)
         call rprintfi1(' '//char(0), i1)
      end subroutine 


      subroutine rprintid3(msg, i1, d1, d2, d3)
      character (len=*) msg
      double precision d1, d2, d3
      integer i1

        call rprintd3(msg//char(0), d1, d2, d3)
        call rprinti1(" "//char(0), i1)
      end subroutine 

C
      subroutine rprintd1(msg, d1)
      character (len=*) msg
      double precision d1
        call rprintfd1(msg//char(0), d1)
      end subroutine 
C
      subroutine rprintd2(msg, d1, d2)
      character (len=*) msg
      double precision d1, d2
        call rprintfd2(msg//char(0), d1, d2)
      end subroutine 
C
      subroutine rprintd3(msg, d1, d2, d3)
      character (len=*) msg
      double precision d1, d2, d3
        call rprintfd3(msg//char(0), d1, d2, d3)
      end subroutine 

C printing with one integer
C printing with integers
      subroutine rprinti1(msg, i1)
      character (len=*) msg
      integer i1
        call rprintfi1(msg//char(0), i1)
      end subroutine 
C
      subroutine rprinti2(msg, i1, i2)
      character (len=*) msg
      INTEGER i1, i2
        call rprintfi2(msg//char(0), i1, i2)
      end subroutine 
C
      subroutine rprinti3(msg, i1, i2, i3)
      character (len=*) msg
      INTEGER i1, i2, i3
        call rprintfi3(msg//char(0), i1, i2, i3)
      end subroutine

C  one logical, one integer
      subroutine rprintli(msg, l1, i1)
      character (len=*) msg
      character (len=8) logi
      logical l1
      integer i1

       if (l1) then
         logi ='  TRUE '
       else
         logi ='  FALSE'
       endif
      
       call rprint(msg//logi//char(0))
       call rprintfi1(msg//char(0), i1)

      end subroutine 
