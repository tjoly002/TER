program main
    use mod_precision
    use mod_maillage
    use mod_sortie

    implicit none
    integer :: CFL,i,k,compteur,nplot
    character(len=40) :: nom_maillage, cl_bord
    integer :: nb_mailles, nb_aretes,nb_arete_bord
    real(pr) :: Delta_t,S,Min,Flux,centre_maille,t_max, h, T_ext, T_init,T_source,D
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


    print*, "il y a ", nb_mailles, " mailles et " ,nb_aretes, " aretes"

    !Initialisation de Tn et Tnp1/Etape3
    allocate(Tn(1:nb_mailles))
    allocate(Tnp1(1:nb_mailles))
    allocate(Residu(1:nb_aretes),Residu_0(1:nb_aretes))

    do i=1,nb_mailles
        Tn(i)=T_init
    end do



    !Calcul du pas de temps si NH/Etape4
     !Calcul du min
    if(cl_bord=="NH")then
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

      print*, "La veleur max du Delta_t pour respecter la CFL vaut : ", Delta_t

      print*, "Quelle valeur de Deltat_t voulez-vous alors ?"

      read*, Delta_t
      print*, Delta_t
    else 
      print*, "Le temps max vaut ", t_max
      print*, "Quelle valeur de Deltat_t voulez-vous pour FR ?"
      read*, Delta_t
      print*, Delta_t
    end if 
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
    !print*, "Boucle en temps 1"
    call sortie(0,Tn,coord_noeud,noeud_maille)
   do i=1,int(t_max/Delta_t)
        !Calcul du flux sur les aretes de bord
        Residu=0.
        do k=1,nb_aretes
            !print*, "Boucle en temps 2"
            Flux=0
            if (maille_arete(k,2)==0) then

                if (cl_arete(k)==10) then
                  if(cl_bord=="FR") then !Condition de Fourier-Robin
                    Flux=-h*D/d_arete(k)*(T_ext-Tn(maille_arete(k,1)))/(h+D/d_arete(k))
                    !print*, "Boucle en temps 3"
                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux
                    end if

                  else if(cl_bord=="NH") then ! Condition de Neumann homogène
                    Flux=0
                    !print*, "Boucle en temps 4"
                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux
                    end if

                  end if

                end if
                if (cl_arete(k)==11) then
                    !print*, "Boucle en temps 5"
                    Flux=-h*D/d_arete(k)*(T_source-Tn(maille_arete(k,1)))/(h+D/d_arete(k))

                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                    else
                      Residu(k)=Flux
                    end if


                else if (cl_arete(k)==12) then
                    !print*, "Boucle en temps 6"
                    Flux=0
                    if (i==1) then
                      Residu_0(k)=Flux
                      Residu(k)=Flux
                      !print*, "Boucle en temps 6.1"
                    else
                      Residu(k)=Flux
                      !print*, "Boucle en temps 6.2"

                    end if

                end if
                !print*, "Boucle en temps 6.3.1", Flux
                Tn(maille_arete(k,1))=Tn(maille_arete(k,1))-Delta_t/aire_maille(maille_arete(k,1))*l_arete(k)*Flux
                !print*, "Boucle en temps 6.3.2"

            else
                !print*, "Boucle en temps 7"
                Flux=-D*(Tn(maille_arete(k,2))-Tn(maille_arete(k,1)))/d_arete(k)
                !print*, "Boucle en temps 7.1"
                Tn(maille_arete(k,1))=Tn(maille_arete(k,1))-Delta_t/aire_maille(maille_arete(k,1))*l_arete(k)*Flux
                !print*, "Boucle en temps 7.2"
                Tn(maille_arete(k,2))=Tn(maille_arete(k,2))+Delta_t/aire_maille(maille_arete(k,2))*l_arete(k)*Flux
                !print*, "Boucle en temps 7.3"

                if (i==1) then
                  Residu_0(k)=Flux
                  Residu(k)=Flux
                else
                  Residu(k)=Flux
                end if

            end if




            Iter(i)=i

        end do

        !print *, "sortie de boucle "



        !print* ,"i=", i
        !print *, "Tn", Tn

        ! Calcul du résidu à chaque itération

        if (i<2000 ) then
         write(2,*), Iter(i),norm2(Residu)/norm2(Residu_0),norm2(Residu)*Delta_t, norm2(Residu)
        !  print*, norm2(Residu)/norm2(Residu_0)
        end if

      ! if(norm2(Residu)/norm2(Residu_0)<0.5) then
      !   print *, 'le nombre d iter max est de ', i
      !   print *, "DeltaT*Residu= ", Delta_t*norm2(Residu)
        
      ! end if



      nplot=int(t_max/Delta_t/100._pr)
      !print*, nplot, t_max
      if (mod(i,nplot)==0) then
        call sortie(i,Tn,coord_noeud,noeud_maille)
        !print*,"Boucle sortie"
       end if
           
    end do
   

    

    !print *, "Tf", Tn

    
  
    deallocate(Tn)
    deallocate(Tnp1)
    deallocate(arete_bord)
    deallocate(Residu)
    deallocate(Iter)
    close(2)

  
end program main
