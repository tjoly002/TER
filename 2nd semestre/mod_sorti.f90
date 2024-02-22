module mod_sortie
    use mod_precision


  implicit none



contains

  subroutine sortie(iter,T,dx,dy)

    ! entré
    integer, intent(in)                :: iter
    real(pr), dimension(:), intent(in) :: T
    real(pr), intent(in)               :: dx,dy

    ! locales
    character(len=30) :: it_ch
    integer           :: i,j, nb_noeuds, nb_noeuds_x, nb_noeuds_y, nb_mailles

    write(it_ch,*) iter
    open(unit=20,file='sortie_'//trim(adjustl(it_ch))//'.vtk')
    !!! A REFLECHIR SURTOUT NIVEAu NX NY AVEC MATRICE AUTRE PART
    nb_noeuds_y = 35/dy+1
    nb_noeuds_x = 0.5/dx+1
    nb_noeuds   = nb_noeuds_x*nb_noeuds_y
    nb_mailles = size(T,1)

    write(20,'(1A26)') '# vtk DataFile Version 2.0'
    write(20,*) 'temperature'
    write(20,*) 'ASCII'
    write(20,*) 'DATASET UNSTRUCTURED_GRID'

    write(20,*) 'POINTS', nb_noeuds,' double'
    ! On remplit ici par ligne
    do i=1,nb_noeuds_y
      do j = 1, nb_noeuds_x
        write(20,*) (i-1)*dy,(j-1)*dx,0.0_pr
      end do
    end do

    write(20,*) 'CELLS ',nb_mailles,nb_mailles*5
    do j=1,nb_noeuds_y-1
      do i = 1,nb_noeuds_x-1
       write(20,*) 4, (j-1)*nb_noeuds_x+i-1, (j-1)*nb_noeuds_x+i, nb_noeuds_x*j+i,nb_noeuds_x*j+i-1
      end do
    end do

    write(20,*) 'CELL_TYPES ', nb_mailles
    do i=1,nb_mailles
       write(20,*) 9
    end do

    write(20,*) 'CELL_DATA',nb_mailles
    write(20,*) 'SCALARS T double'
    write(20,*) 'LOOKUP_TABLE default'
    do i=1,nb_mailles
       write(20,*) T(i)
    end do

    close(20)

  end subroutine sortie

end module mod_sortie
