program main
    use,intrinsic :: iso_fortran_env
    use utils
    implicit none
    
    REAL(real64),PARAMETER :: PI = acos(-1.0_real64)
    INTEGER,PARAMETER :: NUM_TIME_STEP = 5
    INTEGER,PARAMETER :: NUM_POINT = 252, NUM_ELEMENT = 500
    REAL(real64),PARAMETER :: ALMOST0 = 1.0d-8
    REAL(real64),PARAMETER :: WAVE_VELOCITY = 340.0d0
    REAL(real64),PARAMETER :: WAVE_LENGTH = 0.34d0
    REAL(real64),PARAMETER :: TIME_INCREMENT = 0.05d0
    REAL(real64) :: wvwv = WAVE_VELOCITY**2, s_layer, d_layer
    REAL(real64) :: tc

    REAL(real64),ALLOCATABLE :: points(:,:)  ! 点の3軸座標の配列
    INTEGER,ALLOCATABLE :: elements(:,:)  ! 三角形の点番号3つの配列
    INTEGER :: point_num, elem_num, elem_num2, axis_num, elem_point_num, time_step  ! 繰り返し用変数
    REAL(real64) :: element(3,3)  ! 三角形一つの座標を保持する
    REAL(real64) :: center(3)  ! 三角形の重心
    REAL(real64) :: xce(6), yce(3), zce
    REAL(real64),ALLOCATABLE :: yt(:,:,:), yn(:,:,:), yh(:,:)  ! 要素から作る正規直交基底
    REAL(real64),ALLOCATABLE :: s_matrix(:,:,:), d_matrix(:,:,:), c_matrix(:,:)
    REAL(real64),ALLOCATABLE :: elem_u(:,:)
    REAL(real64),ALLOCATABLE :: b_vector(:,:)

    REAL(real64) :: direction(3) = [0.0d0, 0.0d0, 1.0d0]  ! 波の進行方向（単位ベクトル）
    REAL(real64) :: temp, temp2, time
    INTEGER :: boundary_condition(NUM_ELEMENT)  ! 境界条件を要素ごとに設定する配列
    INTEGER :: ipiv(NUM_ELEMENT)  ! lapackで解く際にピボットを保存する配列

    ALLOCATE(points(3,NUM_POINT))
    ALLOCATE(elements(3,NUM_ELEMENT))
    ALLOCATE(yt(3,3,NUM_ELEMENT))
    ALLOCATE(yn, source=yt)
    ALLOCATE(yh(3,NUM_ELEMENT))
    ALLOCATE(elem_u(NUM_ELEMENT, NUM_TIME_STEP), source=0.0d0)
    ALLOCATE(b_vector(NUM_ELEMENT, NUM_TIME_STEP), source=0.0d0)
    ALLOCATE(s_matrix(NUM_ELEMENT, NUM_ELEMENT, 0:NUM_TIME_STEP), source=0.0d0)
    ALLOCATE(d_matrix(NUM_ELEMENT, NUM_ELEMENT, 0:NUM_TIME_STEP), source=0.0d0)
    ALLOCATE(c_matrix(NUM_ELEMENT, NUM_ELEMENT), source=0.0d0)

    ! read mesh data
    open(20, file='mesh/pmesh', status='old')
    open(21, file='mesh/nmesh', status='old')
    do point_num = 1, NUM_POINT
        read(20,*) points(1, point_num), points(2, point_num), points(3, point_num)
    enddo
    do elem_num = 1, NUM_ELEMENT
        read(21,*) elements(1, elem_num), elements(3, elem_num), elements(2, elem_num)
        boundary_condition(elem_num) = -1  ! 1: Dirichlet, -1: Neumann
    end do

    ! write mesh by gnuplot
    call plot(NUM_ELEMENT, points, elements)

    ! 名前付き構文と継続行構文を用いている
    make_boundary_condition :&
    do elem_num = 1, NUM_ELEMENT
        do elem_point_num = 1, 3
            do axis_num = 1, 3
                element(axis_num, elem_point_num) = points(axis_num,elements(elem_point_num, elem_num))
            end do
        end do
        call calc_tnh(element, yt(1,1,elem_num), yn(1,1,elem_num), yh(1,elem_num))

        ! 三角形の重心を計算
        do axis_num = 1, 3
            center(axis_num) = 0.0d0
            do elem_point_num = 1, 3
                center(axis_num) = center(axis_num) + points(axis_num, elements(elem_point_num, elem_num))
            end do
            center(axis_num) = center(axis_num)/3.0d0
        end do

        temp = dot_product(direction, center)/WAVE_VELOCITY
        if (boundary_condition(elem_num) == 1) then  ! Dirichlet
            do time_step = 1, NUM_TIME_STEP
                time = TIME_INCREMENT*time_step
                elem_u(elem_num, time_step)=0.0d0
                if (time > temp) then  ! 波が到達しているなら
                    elem_u(elem_num, time_step) = 1.0d0 - cos(2.0d0*PI*(time - temp)/WAVE_LENGTH)
                endif
            enddo
        else if (boundary_condition(elem_num) == -1) then  ! Neumann
            temp2 = dot_product(direction, yh(:,elem_num))
            temp2 = 2.0d0*pi*temp2/WAVE_VELOCITY/WAVE_LENGTH
            do time_step = 1, NUM_TIME_STEP
                time = TIME_INCREMENT*time_step
                elem_u(elem_num, time_step)=0.0d0
                if (time > temp) then
                    elem_u(elem_num, time_step) = -temp2*sin(2.0d0*PI*(time - temp)/WAVE_LENGTH)
                endif
            enddo
        else
            print *,'wrong b.c.'
            stop
        endif
    end do&
    make_boundary_condition

    ! 先にメモリ空間を0で埋め確保しておく
    ! 割り付け時に0埋めしておいたのでここでしなくても問題ないはず
    ! do time_step = 1, NUM_TIME_STEP
    !     do elem_num = 1, NUM_ELEMENT
    !         b_vector(elem_num, time_step) = 0.0d0;
    !     end do
    ! end do
    ! do time_step = 1, NUM_TIME_STEP
    !     do elem_num = 1, NUM_ELEMENT
    !         do elem_num2 = 1, NUM_ELEMENT
    !             s_matrix(elem_num, elem_num2, time_step) = 0.0d0
    !             d_matrix(elem_num, elem_num2, time_step) = 0.0d0
    !         end do
    !     end do
    ! end do

    !!!! main !!!!
    time_loop :&
    do time_step = 1, NUM_TIME_STEP
        print *,'time_step=',time_step
        time = TIME_INCREMENT*time_step
        tc = time*WAVE_VELOCITY  ! もとは**になっているが*の間違い...?

        element_loop :&
        do elem_num = 1, NUM_ELEMENT
            do elem_point_num = 1, 3
                do axis_num = 1, 3
                    element(axis_num, elem_point_num) = points(axis_num,elements(elem_point_num, elem_num))
                end do
            end do

            point_loop :&
            do point_num = 1, NUM_ELEMENT
                do axis_num = 1, 3
                    center(axis_num) = 0.0d0
                    do elem_point_num = 1, 3
                        center(axis_num) = center(axis_num) + points(axis_num, elements(elem_point_num, point_num))
                    end do
                    center(axis_num) = center(axis_num)/3.0d0
                end do
                CALL calc_xyz(center, element, yt(:,:,elem_num), yn(:,:,elem_num), yh(:,elem_num)&
                                ,xce, yce, zce, ALMOST0)
                
                ! 行列の係数を計算
                s_layer = 0.0d0
                d_layer = 0.0d0
                if ( tc > abs(zce) ) then  ! t*c > |z|
                    do axis_num = 1, 3
                        if (abs(yce(axis_num)) > ALMOST0) then  ! |y| > 0
                            call calc_s_d_layer(tc, xce(2*axis_num-1), yce(axis_num), zce, s_layer, d_layer)
                        endif
                    enddo
                endif
                ! store
                if (boundary_condition(elem_num) == -1) then
                    s_matrix(point_num,elem_num,time_step) = 0.25*s_layer/pi
                    d_matrix(point_num,elem_num,time_step) = 0.25*d_layer/pi
                else if (boundary_condition(elem_num) == 1) then
                    s_matrix(point_num,elem_num,time_step) = -0.25*d_layer/pi
                    d_matrix(point_num,elem_num,time_step) = -0.25*s_layer/pi
                else
                    write(*,*)'wrong b.c.'
                    stop
                endif
            end do point_loop
        end do element_loop

end program main