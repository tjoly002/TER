module mod_sortie

  implicit none

  use mod_precision

contains

  subroutine sortie(iter, T, dx)

    ! entré
    integer, intent(in)                :: iter
    real(pr), dimension(:), intent(in) :: T
    real(pr), intent(in)               :: dx,dy

    ! locales
    character(len=30) :: it_ch
    integer           :: i,j, nb_noeuds, nb_noeuds_x, nb_noeuds_y, nb_mailles

    open(unit=20,file='sortie_'//trim(adjustl(it_ch))//'.vtk')
    nb_noeuds_y = 35._pr%dy
    nb_noeuds_x = 0.5_pr%dx
    nb_noeuds   = (nb_noeuds_x+1)*((nb_noeuds_y+1))
    nb_mailles = size(T,1)

    write(20,'(1A26)') '# vtk DataFile Version 2.0'
    write(20,*) 'temperature'
    write(20,*) 'ASCII'
    write(20,*) 'DATASET UNSTRUCTURED_GRID'

    write(20,*) 'POINTS',nb_noeuds,' double'
    ! On remplit ici par ligne 
    do i=1,nb_noeuds_y
      do j = 1, nb_noeuds_x
        write(20,*) j*dx,i*dy,0.0_pr
      end do
    end do



  end subroutine

end module mod_sortie
