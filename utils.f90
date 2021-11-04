module utils
    use,intrinsic :: iso_fortran_env
    implicit none
contains
    subroutine space_plot(NUM_ELEMENT, points, elements)
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
    end subroutine space_plot

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

        ! calc h
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

    subroutine calc_xyz(point, elem, yt, yn, yh, x, y, z, ALMOST0)
        ! xの1,3,5と2,4,6が各点の座標セットになっている.
        REAL(real64),INTENT(IN) :: point(3), elem(3,3), yt(3,3), yh(3), yn(3,3), ALMOST0
        REAL(real64),INTENT(INOUT) :: x(6), y(3), z
        INTEGER :: axis_num, jm1, k
        REAL(real64) :: temp(3)

        z = dot_product(point(:)-elem(:,1), yh(:))
        if (abs(z) < ALMOST0) z = abs(z)  ! 数値誤差で表裏が入れ替わるのを防ぎ, 小さいときは正であることを保証する
        do axis_num = 1, 3
            ! note: 2*jm1は6,2,4
            jm1=axis_num - 1 + int((4-axis_num)/3)*3
            x(2*axis_num-1) = 0.0d0
            x(2*jm1) = 0.0d0
            y(axis_num) = dot_product(point(:)-elem(:,axis_num), yn(:,axis_num))
            do k = 1, 3
                temp(k)=point(k)-elem(k,axis_num)
                x(2*axis_num-1)=x(2*axis_num-1)+temp(k)*yt(k,axis_num)
                x(2*jm1)=x(2*jm1)+temp(k)*yt(k,jm1)
            end do
        end do
    end subroutine calc_xyz

    subroutine calc_sd_tconst(pm, x, y, z, s_layer,d_layer)
        REAL(real64),INTENT(IN) :: pm, x, y, z
        REAL(real64),INTENT(INOUT) :: s_layer, d_layer
        REAL(real64) :: temp1, temp2, r, ay
        r = sqrt(x**2+y**2+z**2)
        ay = abs(y)

        temp1 = sign(1.0d0, y)*datan2(x*z, ay*r)
        temp2 = sign(1.0d0, y)*datan2(-x, ay)
        s_layer = s_layer + pm*(-z*temp1-y*log(x+r)-abs(z)*temp2)
        d_layer = d_layer + pm*(temp1+sign(1.0d0, z)*temp2)
    
    end subroutine calc_sd_tconst

    subroutine calc_s_d_layer(tc,x,y,z,s_layer,d_layer)
        REAL(real64),INTENT(IN) :: tc, x(2), y, z
        REAL(real64),INTENT(INOUT) :: s_layer, d_layer
        REAL(real64) :: yyzz, tctc, theta, ay, rr, pm, xc
        INTEGER :: i

        yyzz = y**2 + z**2
        tctc = tc**2
        ay = abs(y)
        if (tctc < yyzz) then
            theta = sign(1.0d0, y) * (atan2(-x(2),ay) - atan2(-x(1),ay))
            s_layer = s_layer + (tc-abs(z))*theta
            d_layer = d_layer + sign(1.0d0, z)*theta
            return
        endif

        do i = 1, 2
            rr = x(i)**2 + yyzz
            pm = sign(1.0d0, i - 1.5d0)
            if (tctc < rr) then
                xc = sign(1.0d0, x(i))*sqrt(tctc-yyzz)
                theta = sign(1.0d0, y)*(atan2(-x(i), ay) - atan2(-xc, ay))
                s_layer = s_layer + pm*(tc-abs(z))*theta
                d_layer = d_layer + pm*sign(1.0d0,z)*theta
            else
                xc=x(i)
            endif
            call calc_sd_tconst(pm,xc,y,z,s_layer,d_layer)
        end do
    end subroutine calc_s_d_layer
end module utils
