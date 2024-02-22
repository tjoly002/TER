program main
    use mod_precision
    use mod_sortie

    integer :: iter=1
    real(pr) :: dx=0.1_pr, dy=0.5_pr
    real(pr), dimension(4) :: T

    T=(/0.0,1.0,2.0,3.0/)

    call sortie(iter,T,dx,dy)
end program
