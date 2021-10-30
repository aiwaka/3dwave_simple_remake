module utils
    use,intrinsic :: iso_fortran_env
    implicit none
contains
    subroutine plot(NUM_ELEMENT, points, elements)
        INTEGER,INTENT(IN) :: NUM_ELEMENT, elements(:,:)
        REAL(real64),INTENT(IN) :: points(:,:)
        
        INTEGER :: elem_num
        INTEGER :: j,k

        open(20,file='mesh_gnup.dat')
        ! open(21,file='nor_gnup')
        do elem_num = 1, NUM_ELEMENT
            write(20,*)'#',elem_num,(elements(j,elem_num),j=1,3)
            do j=1,3
            write(20,*)(points(k,elements(j,elem_num)),k=1,3)
            ! do k=1,3
            !     element(k,j)=points(k,elements(j,elem_num))
            ! enddo
            enddo
            write(20,*)(points(k,elements(1,elem_num)),k=1,3)
            write(20,*)
            write(20,*)

            ! call ele_normal(element,vn)
            ! write(21,*)'#',elem_num
            ! write(21,*)(element(k,1)+0.1*vn(k),k=1,3)
            ! write(21,*)(element(k,1),k=1,3)
            ! write(21,*)
            ! write(21,*)
        enddo
        ! close(21)
        close(20)
    end subroutine plot

    subroutine calc_tnh(element, yt, yn, yh)
        REAL(real64),INTENT(IN) :: element(3,3)
        REAL(real64),INTENT(INOUT) :: yt(3,3), yn(3,3), yh(3)
        INTEGER :: i,j,ip1,ip2,jp1
        REAL(real64) :: et, eh, en

        ! calc t
        do j=1,3
            jp1=mod(j,3)+1
            do i=1,3
                yt(i,j)=element(i,jp1)-element(i,j)
            enddo
            et=sqrt(yt(1,j)**2+yt(2,j)**2+yt(3,j)**2)
            do i=1,3
                yt(i,j)=yt(i,j)/et
            enddo
        enddo

        ! calc n
        do i=1,3
            ip1=mod(i,3)+1
            ip2=mod(i+1,3)+1
            yh(i)=-yt(ip1,1)*yt(ip2,3)+yt(ip2,1)*yt(ip1,3)
        enddo
        eh=sqrt(yh(1)**2+yh(2)**2+yh(3)**2)
        do i=1,3
            yh(i)=yh(i)/eh
        enddo
        !calc n
        do j=1,3
            do i=1,3
                ip1=mod(i,3)+1
                ip2=mod(i+1,3)+1
                yn(i,j)=yh(ip1)*yt(ip2,j)-yh(ip2)*yt(ip1,j)
            enddo
            en=sqrt(yn(1,j)**2+yn(2,j)**2+yn(3,j)**2)
            do i=1,3
                yn(i,j)=yn(i,j)/en
            enddo
        enddo        
    end subroutine calc_tnh
end module utils
