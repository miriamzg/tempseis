      integer function lofw(word)

         character*(*) word

         lw=len(word)
         k=0
         do i=1,lw
            if (word(i:i).eq.' ') go to 99
            k=k+1
         end do
  99     lofw=k

         return
      end
