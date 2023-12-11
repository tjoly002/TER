program main
    use mod_precision
    use mod_maillage
    use mod_sortie

    implicit none
    integer :: T_init,T_source,D,CFL,i,k,compteur,nplot
    character(len=40) :: nom_maillage, cl_bord
    integer:: nb_mailles, nb_aretes,nb_arete_bord
    real(pr) :: Delta_t,S,Min,Flux,centre_maille,t_max, h, T_ext
    real(pr), dimension(:), allocatable :: aire_maille, d_arete, l_arete,Tn,Tnp1,Erreur, Residu,Residu_0, Iter
    real(pr), dimension(:,:), allocatable :: coord_noeud, milieu_arete
    integer, dimension(:,:), allocatable :: arete_maille, noeud_maille, maille_arete
    integer, dimension(:), allocatable :: cl_arete,arete_bord






    !Données/Etape1
    open(1,file='donnees')
    read(1,*) T_init
    read(1,*) T_source
    read(1,*) T_ext
    read(1,*) t_max
    read(1,*) cl_bord
    read(1,*) D
    read(1,*) h
    read(1,*) CFL
    read(1,*) nom_maillage


    close(1)
    !print *, T_init,T_G,T_D,Phi_b,t_max,D,CFL,nom_maillage

    open(2,file="Residu.dat")
    !Appel du maillage/Etape2
    call maillage(nom_maillage,nb_mailles, nb_aretes                        &
    &        , coord_noeud, noeud_maille, aire_maille, l_arete, d_arete           &
    &        , milieu_arete, arete_maille, maille_arete, cl_arete)




    !Initialisation de Tn et Tnp1/Etape3
    allocate(Tn(1:nb_mailles))
    allocate(Tnp1(1:nb_mailles))
    allocate(Residu(1:nb_aretes),Residu_0(1:nb_aretes))

    do i=1,nb_mailles
        Tn(i)=T_init
    end do



    !Calcul du pas de temps /Etape4
     !Calcul du min
    min=10e6
    do i=1,nb_mailles

        !Calcul de la somme pour Delta_t

        S=0
        do k=1,3

            S=S+l_arete(arete_maille(i,k))*D/d_arete(arete_maille(i,k))
        end do

        if (aire_maille(i)/S<min) then
            min=aire_maille(i)/S

        end if

    end do

    Delta_t=min



    allocate(Iter(1:int(t_max/Delta_t)))

    !Boucle en temps/Etape5

    !Calcul nombre d'aretes de bords
    nb_arete_bord=0
    do k=1,nb_aretes
        if (maille_arete(k,2)==0) then
            nb_arete_bord=nb_arete_bord+1
        end if
    end do

    !Allocate le vecteur contenant les numeros des aretes de bord
    allocate(arete_bord(nb_arete_bord))

    compteur=1 !indice du vecteur contenant les numeros des aretes de bord

    do k=1,nb_aretes
        if (maille_arete(k,2)==0) then
            arete_bord(compteur)=k
            compteur=compteur+1
        end if
    end do


    !Boucle en temps
  !  print*, "Boucle en temps 1"
    call sortie(0,Tn,coord_noeud,noeud_maille)
   do i=1,int(t_max/Delta_t)
        !Calcul du flux sur les aretes de bord
        Residu=0.
        do k=1,nb_aretes

            Flux=0
            if (maille_arete(k,2)==0) then

                if (cl_arete(k)==10) then
                  if(cl_bord=="FR") then !Condition de Fourier-Robin
                    Flux=h*(Tn(maille_arete(k,1))-T_ext)

                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux
                    end if

                  else if(cl_bord=="NH") then ! Condition de Neumann homogène
                    Flux=0

                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux
                    end if

                  end if

                end if
                if (cl_arete(k)==11) then

                    Flux=-1._pr*D*(T_source-Tn(maille_arete(k,1)))/d_arete(k)

                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux
                    end if


                else if (cl_arete(k)==12) then

                    Flux=0
                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux

                    end if

                end if

                Tn(maille_arete(k,1))=Tn(maille_arete(k,1))-Delta_t/aire_maille(maille_arete(k,1))*l_arete(k)*Flux
            else

                Flux=-D*(Tn(maille_arete(k,2))-Tn(maille_arete(k,1)))/d_arete(k)
                Tn(maille_arete(k,1))=Tn(maille_arete(k,1))-Delta_t/aire_maille(maille_arete(k,1))*l_arete(k)*Flux
                Tn(maille_arete(k,2))=Tn(maille_arete(k,2))+Delta_t/aire_maille(maille_arete(k,2))*l_arete(k)*Flux

                if (i==1) then
                  Residu_0(k)=Flux
                  Residu(k)=Flux
                else
                  Residu(k)=Flux
                end if

            end if




            Iter(i)=i

        end do
        ! print* ,"i=", i
        ! print *, "Tn", Tn

        ! Calcul du résidu à chaque itération
      if (i>100 .and. i<200) then
      write(2,*), Iter(i),norm2(Residu)/norm2(Residu_0)
      end if

      ! if(norm2(Residu)/norm2(Residu_0)<0.123) then
      !   print *, 'le nombre d iter max est de ', i
      !   print *, "DeltaT*Residu= ", Delta_t*norm2(Residu)
      !   exit
      ! end if




        nplot=int(t_max/Delta_t/100)
        if (mod(i,nplot)==0) then
            call sortie(i,Tn,coord_noeud,noeud_maille)
        end if
    end do
    !print*, norm2(Residu_0)

    !print *, "Tf", Tn

    !Calcul de l'erreur entre notre résultat et la fonction exacte
    ! allocate(Erreur(nb_mailles))
    ! do i=1,nb_mailles
    !     centre_maille=(coord_noeud(noeud_maille(i,1),1)+coord_noeud(noeud_maille(i,2),1)+coord_noeud(noeud_maille(i,3),1))/3
    !     !print *, centre_maille, coord_noeud(noeud_maille(i,1),1), coord_noeud(noeud_maille(i,3),1)
    !     Erreur(i)=abs(Tn(i)-fonction_exacte(1._pr*t_max,centre_maille))
    ! end do
    ! print*, Normeinf(Erreur)

  !   deallocate(Erreur)
    deallocate(Tn)
    deallocate(Tnp1)
    deallocate(arete_bord)
    deallocate(Residu)
    deallocate(Iter)
    close(2)

    !
    ! contains
    ! !fonction exacte
    ! function fonction_exacte(t,r) result(y)
    !     real(pr), intent(in) :: t, r
    !     real(pr) :: y
    !
    !     integer :: n
    !     real(pr), dimension(1:20) :: z, c
    !
    !     !--- zeros de J0
    !     z = (/ 2.404825558, 5.520078110, 8.653727913, 11.79153444, 14.93091771 &
    !         & , 18.07106397, 21.21163663, 24.35247153, 27.49347913, 30.63460647 &
    !         & , 33.77582021, 36.91709835, 40.05842576, 43.19979171, 46.34118837 &
    !         & , 49.48260990, 52.62405184, 55.76551076, 58.90698393, 62.04846919/)
    !
    !     !--- coeff cn
    !     do n=1,20
    !         c(n) = 2 * (T_init - T_source) / ( z(n) * Bessel_JN(1,z(n)) )
    !     end do
    !
    !     !--- T exacte
    !     y  = T_source
    !     do n=1,20
    !         y = y + c(n) * exp(-z(n)**2 * D/**2 * t) * Bessel_JN(0,r*z(n)/ )
    !     end do
    !
    ! end function fonction_exacte
    !
    ! !foncion calculant la norme 2
    ! function norme2(Vect) result(res)
    !     real(pr),dimension(:),intent(in) :: Vect
    !     real(pr) :: res
    !     integer :: i
    !
    !     res=0
    !     do i=1,size(Vect)
    !         res=res+Vect(i)**2
    !     end do
    !     res=sqrt(res)
    ! end function Norme2
    !
    ! !foncion calculant la norme inf
    ! function normeinf(Vect) result(res)
    !     real(pr),dimension(:),intent(in) :: Vect
    !     real(pr) :: res
    !     integer :: i
    !
    !     res=0
    !     do i=1,size(Vect)
    !         if (res<Vect(i)) then
    !             res=Vect(i)
    !         end if
    !     end do
    ! end function Normeinf
end program main
