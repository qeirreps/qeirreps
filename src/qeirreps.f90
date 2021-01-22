program Main
  implicit none
  integer :: cmdcount, Nk_irr, NBAND, ncomp, nsymq, nnp, &
       ik, ib, ic, ig, iop, i, j, &
       maxNGp, maxNGI, tmpkg(3), iinv, filling
  real(8) :: Ecut_for_psi, FermiEnergy, Etot, a1(3), a2(3), a3(3)
  logical :: flag,flag_fill,flag_z4,flag_ctm
  character(150) :: data_dir,data_fill,data_option
  integer, allocatable ::  NGI(:), KGI(:,:,:), rg(:,:,:), pg(:,:),product_list(:,:),&
       NG_for_psi(:), tfkg(:,:,:), nband_occ(:), &
       degeneracy(:,:), energy_level(:),Gshift(:,:,:)
  real(8), allocatable :: SKI(:,:), E_EIGI(:,:), tr(:,:,:),rotA(:,:,:)
  complex(8), allocatable :: C0(:,:,:,:),repmat(:,:,:),chrct(:,:,:), &
       sr(:,:,:),phase(:,:,:),kphase(:,:),facsys_spin(:,:,:)
  logical, allocatable :: gk(:,:)

  maxNGp = 0
  maxNGI = 0


!============================================

  ! get data

!============================================

  ! get the name of the directory, filling, and the option
  flag_fill=.false.
  flag_z4=.false.
  flag_ctm=.false.

  cmdcount=command_argument_count()
  if(cmdcount == 0)then
     stop "please specify a directory"
  else if(cmdcount == 1)then
     flag_fill=.false.
     flag_z4=.false.
     flag_ctm=.false.
     call get_command_argument(1,data_dir)
     print*,"not specify the filling"
     print*,"reading:"//trim(data_dir)
  else if(cmdcount == 2)then
     flag_fill=.true.
     flag_z4=.true.
     flag_ctm=.false.
     call get_command_argument(1,data_dir)
     call get_command_argument(2,data_fill)
     read(data_fill,*)filling
     print*,"reading:"//trim(data_dir)
     print*,"filling:",filling
     print*,"calculating the Z4 index"
  else if(cmdcount == 3)then
     flag_fill=.true.
     call get_command_argument(1,data_dir)
     call get_command_argument(2,data_fill)
     call get_command_argument(3,data_option)
     read(data_fill,*)filling
     print*,"reading:"//trim(data_dir)
     print*,"filling:",filling
     if(trim(data_option)=="z4")then
        flag_z4=.true.
        flag_ctm=.false.
        print*,"calculating the Z4 index"
     else if(trim(data_option)=="ctm")then
        flag_z4=.false.
        flag_ctm=.true.
        print*,"exporting an input file for CheckTopologicalMat"
     else
        stop "Option error: Please specify 'z4' or 'ctm.'"
     end if
        
  else
     stop "Reading Error: Please specify the directory name, filling, and the option"
  end if
  

  ! reading data from the directory "dir-wfn" which contains
  ! the output files of qe2respack.py
  OPEN(117,FILE=trim(data_dir)//"/dir-wfn/dat.bandcalc")
  read(117,*) Ecut_for_psi   ! cutoff energy for wavefunction (Ry)
  read(117,*) FermiEnergy    ! Fermi energy (au) real(8)
  read(117,*) Etot           ! total energy (au) real(8)
  CLOSE(117)


  OPEN(101,FILE=trim(data_dir)//"/dir-wfn/dat.sample-k")
  read(101,*) Nk_irr   ! number of sample k-points
  allocate(SKI(3,Nk_irr))
  do ik=1,Nk_irr
     read(101,*) (SKI(i,ik),i=1,3) ! sample k-points
  enddo
  CLOSE(101)


  OPEN(111,FILE=trim(data_dir)//"/dir-wfn/dat.eigenvalue")
  rewind(111)
  read(111,*) NBAND   ! number of bands
  allocate(E_EIGI(NBAND,Nk_irr))
  do ik=1,Nk_irr
     do ib=1,NBAND
        read(111,*) E_EIGI(ib,ik) ! band energy (au)
     end do
  end do
  CLOSE(111)


  OPEN(132,FILE=trim(data_dir)//"/dir-wfn/dat.nkm")
  rewind(132)
  allocate(NGI(Nk_irr)) ! number of expansion reciprocal lattice vectors
  do ik=1,Nk_irr
     read(132,*) NGI(ik) 
  end do
  CLOSE(132)


  OPEN(104,FILE=trim(data_dir)//"/dir-wfn/dat.kg")
  allocate(NG_for_psi(Nk_irr))
  do ik=1,Nk_irr
     read(104,*) NG_for_psi(ik) ! number of expansion reciprocal lattice vectors
     do ig =1,NG_for_psi(ik)
        read(104,'()')
     end do
  end do

  rewind(104)
  maxNGp = maxval(NG_for_psi)
  allocate(KGI(3,maxNGp,Nk_irr))
  
  do ik=1,Nk_irr
     read(104,'()')
     do ig=1,NG_for_psi(ik)
        read(104,*) (KGI(i,ig,ik),i=1,3) ! expansion reciprocal lattice vectors
     end do
  end do
  CLOSE(104)


  OPEN(102,FILE=trim(data_dir)//"/dir-wfn/dat.wfn",FORM='unformatted')
  rewind(102)
  read(102) ncomp ! number of components of wavefunction
                  ! ncomp = 2 for calculation with spin-orbit coupling
  
  print*,"ncomp=",ncomp
  print*,""

  maxNGI = maxval(NGI)
  allocate(C0(ncomp,maxNGI,NBAND,Nk_irr))
  C0=0.0d0

  do ik=1,Nk_irr
     do ib =1,NBAND
        read(102) ((C0(ic,i,ib,ik),i=1,NGI(ik)),ic=1,ncomp) ! wavefunction
     end do
  end do

  CLOSE(102)


  OPEN(105,FILE=trim(data_dir)//"/dir-wfn/dat.lattice")
  read(105,*) a1(1), a1(2), a1(3) ! lattice vector 1
  read(105,*) a2(1), a2(2), a2(3) ! lattice vector 2
  read(105,*) a3(1), a3(2), a3(3) ! lattice vector 3
  CLOSE(105)


  OPEN(100,FILE=trim(data_dir)//"/dir-wfn/dat.symmetry")
  read(100,*) nsymq ! number of symmetry
  read(100,*) nnp ! scale of fractional translation
                  ! fractional translation can be written as N/nnp
  allocate(rg(3,3,nsymq))
  allocate(pg(3,nsymq))
  do iop=1,nsymq
     read(100,*) ((rg(i,j,iop),i=1,3),j=1,3) ! point group part of symmetry
                                             ! in the primitive coordinate
                                             ! of the reciprocal space
     read(100,*) (pg(i,iop),i=1,3)           ! translation part of symmetry
                                             ! with the primitive coordinate
  end do
  CLOSE(100)

!============================================

  ! Main

!============================================

  ! output the reciprocal lattice vectors
  call get_reci_vec(a1,a2,a3)
  
  ! output the operation type for each group operator
  call operation_type(nsymq,rg,pg,iinv)

  allocate(tr(3,3,nsymq)) ! point group part of symmetry operators
                          ! in the conventional coordinate
  allocate(rotA(3,3,nsymq)) ! point group part of symmetry operators
                            ! in the primitive coordinate
                            ! of the real space
  
  ! output the point group part of symmetry operators
  ! in the conventional coordinate
  call rg_out(nsymq,rg,a1,a2,a3,tr,rotA,data_dir)

  ! output the translation part of symmetry operators
  ! in the primitive coordinate
  call pg_out(nsymq,nnp,pg,data_dir)

  allocate(sr(2,2,2*nsymq))

  ! get and output the spin rotation part of symmetry operators
  call spin_rot(nsymq,tr,sr,data_dir)

  allocate(nband_occ(Nk_irr))

  ! specify the number of bands for calculation
  call get_occ_bands_full(NBAND,nband_occ)

  allocate(energy_level(Nk_irr))     ! number of energy level for each k
  allocate(degeneracy(NBAND,Nk_irr)) ! number of degeneracy for each level

  ! get how degeneracy occurs among considered bands
  call get_degeneracy_full_levels(Nk_irr,E_EIGI,NBAND,degeneracy,energy_level)


  allocate(gk(Nk_irr,nsymq)) ! whether the symmetry operator is included in k-group
  allocate(Gshift(3,Nk_irr,nsymq)) ! shift of k-points by symmetry operation
  allocate(kphase(Nk_irr,nsymq)) ! the phase factor of representations in form of
                                 ! Exp(-ik t_g)

  ! get k-group G_k and the phase factor of represenatations
  call get_kgroup(Nk_irr,nsymq,SKI,rg,pg,nnp,gk,Gshift,kphase) 


  allocate(tfkg(maxNGp,Nk_irr,nsymq)) ! indices of expansion reciprocal lattice vectors
                                      ! transformed by symmetry operation
  allocate(phase(maxNGp,Nk_irr,nsymq)) ! the phase factor of representation in form of
                                       ! Exp(-iG t_g)

  ! get map of expansion reciprocal lattice vectors by symmetry operation
  ! and phase factor of representations
  call get_transformed_G(Nk_irr,nsymq,maxNGp,NG_for_psi,KGI,gk,Gshift,pg,nnp,tfkg,phase)

  allocate(repmat(NBAND,Nk_irr,nsymq)) ! diagonal part of the representation matrix

  ! get diagonal part of the representation matrix
  call representation_matrix(Nk_irr,nsymq,ncomp,maxNGI,NBAND,nband_occ,gk,NGI,kphase,tfkg,phase,sr,C0,repmat)


  allocate(facsys_spin(Nk_irr,nsymq,nsymq)) ! spin factor systems
  allocate(product_list(nsymq,nsymq)) ! product table of group operators
  
  ! get and output spin factor systems
  call get_spin_factor_system(Nk_irr,nsymq,tr,gk,sr,facsys_spin,product_list)
  call factor_system_out(Nk_irr,nsymq,gk,SKI,facsys_spin,data_dir)

  
  allocate(chrct(Nk_irr,nsymq,maxval(energy_level))) ! characters

  ! get and output characters
  call get_character(Nk_irr,nsymq,NBAND,energy_level,degeneracy,repmat,chrct)
  call character_out(Nk_irr,nsymq,energy_level,SKI,gk,chrct,data_dir)

  ! get and output character tables of irreducible representations for each k
  call character_table_out(Nk_irr,nsymq,ncomp,energy_level,&
       rotA,nnp,pg,SKI,gk,chrct,product_list,facsys_spin,data_dir)

  
  ! get and output z4 index
  ! this process runs only when filling is specified
  if(flag_fill .and. flag_z4) then
     call get_z4(Nk_irr,nsymq,NBAND,iinv,filling,gk,energy_level,degeneracy,chrct,data_dir)
  end if
  
  ! output the file for CheckTopologicalMat
  ! this process runs only when filling is specified
  
  if(flag_fill .and. flag_ctm) then
     call Bilbao_out(Nk_irr,nsymq,ncomp,NBAND,filling,energy_level,degeneracy,&
       E_EIGI,nnp,pg,sr,rotA,SKI,gk,chrct,data_dir)
  end if
  
  print*,"================="
  print*,""
  print*,"job done"

!============================================

  ! End

!============================================

  deallocate(NGI); deallocate(KGI); deallocate(rg); deallocate(pg)
  deallocate(SKI); deallocate(E_EIGI)
  deallocate(C0)
  deallocate(nband_occ)
  deallocate(gk); deallocate(Gshift); deallocate(kphase)
  deallocate(tr); deallocate(sr)
  deallocate(degeneracy); deallocate(energy_level)
  deallocate(facsys_spin)
  deallocate(repmat)
  deallocate(chrct)


contains

!============================================

  ! Subroutines

!============================================

  ! output the reciprocal lattice vectors
  subroutine get_reci_vec(a1,a2,a3)
    implicit none
    real(8),intent(in) :: a1(3),a2(3),a3(3) ! lattice vectors

    real(8),parameter :: pi=acos(-1d0)
    real(8) :: vol,b(3,3)

    ! volume of lattice
    vol =abs(a1(1)*a2(2)*a3(3) - a1(1)*a2(3)*a3(2)&
         & + a1(2)*a2(3)*a3(1) - a1(2)*a2(1)*a3(3)&
         & + a1(3)*a2(1)*a3(2) - a1(3)*a2(2)*a3(1))

    print*,"================="
    print*,""
    
    print*,"lattice vec: bohr"
    print*,a1
    print*,a2
    print*,a3
    print*,""

    ! reciprocal lattice vectors
    b(1,1)=2*pi*(a2(2)*a3(3)-a2(3)*a3(2))/vol
    b(1,2)=2*pi*(a2(3)*a3(1)-a2(1)*a3(3))/vol
    b(1,3)=2*pi*(a2(1)*a3(2)-a2(2)*a3(1))/vol

    b(2,1)=2*pi*(a3(2)*a1(3)-a3(3)*a1(2))/vol
    b(2,2)=2*pi*(a3(3)*a1(1)-a3(1)*a1(3))/vol
    b(2,3)=2*pi*(a3(1)*a1(2)-a3(2)*a1(1))/vol

    b(3,1)=2*pi*(a1(2)*a2(3)-a1(3)*a2(2))/vol
    b(3,2)=2*pi*(a1(3)*a2(1)-a1(1)*a2(3))/vol
    b(3,3)=2*pi*(a1(1)*a2(2)-a1(2)*a2(1))/vol

    print*,"reci_vec: 1/bohr"
    print*,b(1,:)
    print*,b(2,:)
    print*,b(3,:)
    print*,""
    
  end subroutine get_reci_vec
  
  ! output the operation type for each group operator
  subroutine operation_type(nsymq,rg,pg,iinv)
    implicit none
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: rg(3,3,nsymq) ! point group part of symmetry operation
    integer,intent(in) :: pg(3,nsymq) ! translation part of symmetry operation
    
    integer,intent(out) :: iinv ! index of inversion symmetry operator
    integer :: i, j, iop
    real(8) :: err, tmp1(3,3),tmp2(3,3), id(3,3), inv(3,3), det

    err = 0.00001
    id = 0d0
    inv = 0d0
    iinv=0
    
    ! identity operator
    id(1,1)=1d0; id(2,2)=1d0; id(3,3)=1d0

    ! inversion operator
    inv(1,1)=-1d0; inv(2,2)=-1d0; inv(3,3)=-1d0

    print*,"================="
    print*,"symmetry operation type:"

    do iop=1,nsymq

       ! check determinant of the matix
       det=rg(1,1,iop)*rg(2,2,iop)*rg(3,3,iop)&
            &+rg(1,2,iop)*rg(2,3,iop)*rg(3,1,iop)&
            &+rg(1,3,iop)*rg(2,1,iop)*rg(3,2,iop)&
            &-rg(1,1,iop)*rg(2,3,iop)*rg(3,2,iop)&
            &-rg(1,2,iop)*rg(2,1,iop)*rg(3,3,iop)&
            &-rg(1,3,iop)*rg(2,2,iop)*rg(3,1,iop)

       ! check and output the operation type of each symmetry operator

       ! case of identity
       if(maxval(abs(id(:,:)-rg(:,:,iop)))<err)then
          if(maxval(abs(pg(:,iop)))>0)then
             print*,iop,"translation"
          else
             print*,iop,"E"
          end if

       ! case of inversion
       else if(maxval(abs(inv(:,:)-rg(:,:,iop)))<err)then
          iinv=iop
          if(maxval(abs(pg(:,iop)))>0)then
             print*,iop,"I+translation"
          else
             print*,iop,"I"             
          end if
       else if(det>0d0)then
          tmp1(:,:)=rg(:,:,iop)
          do i=1,3
             do j=1,3
                tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
             end do
          end do

          ! case of C2 or 2-fold screw
          if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
             if(maxval(abs(pg(:,iop)))>0)then
                print*,iop,"2-fold screw"
             else
                print*,iop,"C2"
             end if
          else
             tmp1(:,:)=tmp2(:,:)
             do i=1,3
                do j=1,3
                   tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                end do
             end do

             ! case of C3 or 3-fold screw
             if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
                if(maxval(abs(pg(:,iop)))>0)then
                   print*,iop,"3-fold screw"
                else
                   print*,iop,"C3"
                end if
             else
                tmp1(:,:)=tmp2(:,:)
                do i=1,3
                   do j=1,3
                      tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                   end do
                end do

                ! case of C4 or 4-fold screw
                if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
                   if(maxval(abs(pg(:,iop)))>0)then
                      print*,iop,"4-fold screw"
                   else
                      print*,iop,"C4"
                   end if
                else
                   tmp1(:,:)=tmp2(:,:)
                   do i=1,3
                      do j=1,3
                         tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                      end do
                   end do

                   ! case of C6 or 6-fold screw
                   if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
                      if(maxval(abs(pg(:,iop)))>0)then                         
                         print*,iop,"6-fold screw"
                      else
                         print*,iop,"C6"
                      end if
                   end if
                end if
             end if
          end if

       else if(det<0d0)then
          tmp1(:,:)=rg(:,:,iop)
          do i=1,3
             do j=1,3
                tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
             end do
          end do

          ! case of mirror or glide
          if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
             if(maxval(abs(pg(:,iop)))>0)then
                print*,iop,"glide"
             else
                print*,iop,"M"
             end if
          else
             tmp1(:,:)=tmp2(:,:)
             do i=1,3
                do j=1,3
                   tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                end do
             end do
             
             ! case of IC3 or IC3+translation
             if(maxval(abs(inv(:,:)-tmp2(:,:)))<err)then
                if(maxval(abs(pg(:,iop)))>0)then
                   print*,iop,"IC3+translation"
                else
                   print*,iop,"IC3"
                end if
             else
                tmp1(:,:)=tmp2(:,:)
                do i=1,3
                   do j=1,3
                      tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                   end do
                end do

                ! case of IC4 or IC4+translation
                if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
                   if(maxval(abs(pg(:,iop)))>0)then
                      print*,iop,"IC4+translation"
                   else
                      print*,iop,"IC4"
                   end if
                else
                   tmp1(:,:)=tmp2(:,:)
                   do i=1,3
                      do j=1,3
                         tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                      end do
                   end do

                   tmp1(:,:)=tmp2(:,:)
                   do i=1,3
                      do j=1,3
                         tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                      end do
                   end do
                   
                   ! case of IC6 or IC6+translation
                   if(maxval(abs(id(:,:)-tmp2(:,:)))<err)then
                      if(maxval(abs(pg(:,iop)))>0)then
                         print*,iop,"IC6+translation"
                      else
                         print*,iop,"IC6"
                      end if
                   end if
                end if
             end if
          end if

       end if

    end do
    print*,""
    
  end subroutine operation_type

  ! output the point group part of symmetry operators
  ! with the conventional coordinate
  subroutine rg_out(nsymq,rg,a1,a2,a3,tr,rotA,data_dir)
    implicit none
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: rg(3,3,nsymq) ! point group part of symmetry operators
                                        ! in the primitive coordinate
                                        ! of the reciprocal space 
    real(8),intent(in) :: a1(3),a2(3),a3(3) ! lattice vectors
    real(8),intent(out) ::  tr(3,3,nsymq) ! rotation part of symmetry operators
                                          ! with the conventional coordinate
    real(8),intent(out) :: rotA(3,3,nsymq) ! rotation part of symmetry operators
                                           ! with the primitive coorditnate
                                           ! of the real space
    character(150),intent(in) :: data_dir ! output directory

    real(8) :: lat(3,3), latinv(3,3), tmp1(3,3)
    real(8) :: MatA(3,3),MatB(3,3)

    integer :: i, j, iop, info, lda, lwork, n=3
    real(8),allocatable :: a(:,:), work(:)
    integer,allocatable :: ipiv(:)

    lda = n
    lwork = 64*n
    allocate(a(lda,n),work(lwork),ipiv(n))

    lat(:,1)=a1(:)
    lat(:,2)=a2(:)
    lat(:,3)=a3(:)

    a(:,:)=lat(:,:)


    ! get the inversion matrix of lat(3,3)
    call dgetrf(n,n,a,lda,ipiv,info)
    if(info==0)then
       call dgetri(n,a,lda,ipiv,work,lwork,info)
       latinv(:,:)=a(:,:)
    else
       print*,"The factor U is singular"
    end if


    ! get point group reperesentation in the conventional coordinate
    ! by calculating lat^-1 rg lat
    do iop=1,nsymq

       do i=1,3
          do j= 1,3
             tmp1(i,j)=sum(rg(i,:,iop)*lat(j,:))
          end do
       end do

       do i=1,3
          do j= 1,3
             tr(i,j,iop)=sum(latinv(:,i)*tmp1(:,j))
          end do
       end do
       
    end do

    ! output the point group part of symmetry operator in the conventional coordinate
    OPEN(11,FILE=trim(data_dir)//"/output/pg.dat",status="replace")

    do iop=1,nsymq
       write(11,*)  iop
       write(11,*) tr(1,:,iop)
       write(11,*) tr(2,:,iop)
       write(11,*) tr(3,:,iop)
       write(11,*) ""
    end do
    
    CLOSE(11)

    OPEN(12,FILE=trim(data_dir)//"/output/pg_import.txt",status="replace")

    write(12,fmt='(a)',advance='no')"{"
    do iop=1,nsymq
       write(12,'("{{",1f12.8,",",1f12.8,",",1f12.8,"}, &
            & {",1f12.8,",",1f12.8,",",1f12.8,"},{",1f12.8,",",1f12.8,",",1f12.8,"}}")' &
            & ,advance='no') tr(1,:,iop),tr(2,:,iop),tr(3,:,iop)
       if(iop/=nsymq)then
          write(12,fmt='(a)',advance='no')","
       end if
    end do
    write(12,fmt='(a)',advance='no')"}"
    
    CLOSE(12)

    ! get rotation part of symmetry operators with the primitive coorditnate of the real space

    ! transformation matrix between primitive basis of the reciprocal space
    ! and the primitive basis of the real space
    do i = 1,3
       do j = 1,3
          MatA(i,j)=sum(lat(:,i)*lat(:,j))
       end do
    end do

    a(:,:)=MatA(:,:)

    call dgetrf(n,n,a,lda,ipiv,info)

    if(info==0)then
       call dgetri(n,a,lda,ipiv,work,lwork,info)
       MatB(:,:)=a(:,:)

    else
       print*,"Error for transforming reci. lat. to lat."
    end if

    do iop=1,nsymq
       
       tmp1=0d0
       do i=1,3
          do j=1,3
             tmp1(i,j)=sum(rg(i,:,iop)*MatA(:,j))
          end do
       end do
       
       a=0d0
       do i=1,3
          do j=1,3
             rotA(i,j,iop)=sum(MatB(i,:)*tmp1(:,j))
          end do
       end do
    end do

    deallocate(a,work,ipiv)
    
  end subroutine rg_out


  ! output the translation part of symmetry operators
  ! with the primitive coordinate
  subroutine pg_out(nsymq,nnp,pg,data_dir)
    implicit none
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: nnp ! scale of fractional translation
                              ! fractional translation can be written as N/nnp
    integer,intent(in) :: pg(3,nsymq) ! translation part of symmetry operators
    character(150),intent(in) :: data_dir ! output directory
    
    integer :: iop


    ! output the translation part of symmetry operators
    ! with the primitive coordinate
    OPEN(11,FILE=trim(data_dir)//"/output/tg.dat",status="replace")

    do iop=1,nsymq
       write(11,*)iop
       write(11,'(1I3,"/",1I3," ",1I3,"/",1I3," ",1I3,"/",1I3)')pg(1,iop),nnp,pg(2,iop),nnp,pg(3,iop),nnp
       write(11,*)""
    end do

    CLOSE(11)
    
    OPEN(12,FILE=trim(data_dir)//"/output/tg_import.txt",status="replace")

    write(12,fmt='(a)',advance='no')"{"
    do iop=1,nsymq
       write(12,'("{",1I3,"/",1I3,",",1I3,"/",1I3,",",1I3,"/",1I3,"}")' &
            ,advance='no') pg(1,iop),nnp,pg(2,iop),nnp,pg(3,iop),nnp
       if(iop/=nsymq)then
          write(12,fmt='(a)',advance='no')","
       end if
    end do
    write(12,fmt='(a)',advance='no')"}"
    
    CLOSE(12)

  end subroutine pg_out

  ! get and output the spin rotation part of symmetry operators
  subroutine spin_rot(nsymq,tr,sr,data_dir)
    implicit none
    integer,intent(in) :: nsymq ! number of symmetry operation
    real(8),intent(in) ::  tr(3,3,nsymq) ! point group part of symmetry operation
                                         ! in the conventional coordinate
    complex(8),intent(out) ::  sr(2,2,2*nsymq) ! spin rotation matrix
    character(150),intent(in) :: data_dir ! output directory

    complex,parameter :: ii=(0,1)
    integer :: i, j, iop, n=3, lda, lwork, info
    real(8) :: err, id(3,3), inv(3,3), det, vec(3), nvec(3), norm, trace, theta
    real(8),allocatable :: w(:),rwork(:)
    complex(8) :: tmp1(3,3),tmp2(3,3), wq(1)
    complex(8),allocatable :: a(:,:), work(:)

    lda=n
    allocate(a(n,n))
    allocate(w(n));
    allocate(rwork(5*n))


    sr=0d0
    
    err = 0.00001
    id = 0d0
    inv = 0d0

    tmp1=0d0; tmp2=0d0

    ! 3 by 3 identity matrix
    id(1,1)=1d0; id(2,2)=1d0; id(3,3)=1d0

    ! 3 by 3 inversion matrix
    inv(1,1)=-1d0; inv(2,2)=-1d0; inv(3,3)=-1d0

    
    do iop=1,nsymq
       tmp1=0d0; tmp2=0d0
       vec=0d0; nvec=0d0; norm=0d0; trace=0d0; theta=0

       ! get determination of the rotation matrix
       det=tr(1,1,iop)*tr(2,2,iop)*tr(3,3,iop)&
            &+tr(1,2,iop)*tr(2,3,iop)*tr(3,1,iop)&
            &+tr(1,3,iop)*tr(2,1,iop)*tr(3,2,iop)&
            &-tr(1,1,iop)*tr(2,3,iop)*tr(3,2,iop)&
            &-tr(1,2,iop)*tr(2,1,iop)*tr(3,3,iop)&
            &-tr(1,3,iop)*tr(2,2,iop)*tr(3,1,iop)

       ! if the point group part include inversion operation as I R,
       ! map it to the rotation operation as R
       if(abs(det+1)<err)then
          tmp1(:,:)=-tr(:,:,iop)
       else if(abs(det-1)<err)then
          tmp1(:,:)=tr(:,:,iop)
       end if

       do i=1,n
          do j=1,n
             tmp2(i,j)=tmp1(i,j)-tmp1(j,i)
          end do
       end do

       ! get rotation axis
       
       ! the case of identity matrix
       if(maxval(abs(tmp1(:,:)-id(:,:)))<err)then
          sr(:,:,iop)=reshape((/1d0,0d0,0d0,1d0/),(/2,2/))
          cycle
       ! the case of symmetric matrix
       else if(maxval(abs(tmp2(:,:)))<err)then

          a(:,:)=tmp1(:,:)

          lwork = -1
          call zheev("V","U",n,a,lda,w,wq,lwork,rwork,info)
          lwork = Int(wq(1))
          allocate(work(lwork)); work = 0d0
          
          call zheev ("V","U",n,a,lda,w,work,lwork,rwork,info) ! diagonalize 
     
          deallocate(work)

          do i=1,n
             if(abs(w(i)-1)<err)then
                vec(:)=a(:,i)
             end if
          end do
       ! other cases
       else
          vec(1)=tmp2(3,2)
          vec(2)=tmp2(1,3)
          vec(3)=tmp2(2,1)
       end if

       ! normalize the vector of rotation axis
       norm=sqrt(sum(vec(:)**2))
       nvec(:)=vec(:)/norm

       ! get rotation angle
       trace=tmp1(1,1)+tmp1(2,2)+tmp1(3,3)

       if(trace/2d0-5d-1>=-1)then
          if(trace/2d0-5d-1<=1)then
             theta=acos(trace/2d0-5d-1)
          else
             print*,'set theta=pi at',iop
             theta=acos(1d0)
          end if
       else
          print*,'set theta=-pi at',iop
          theta=acos(-1d0)
       end if

       ! get spin rotation matrix from rotation axis and angle
       sr(1,1,iop)=cos(theta/2d0)-ii*nvec(3)*sin(theta/2d0)
       sr(1,2,iop)=-(ii*nvec(1)+nvec(2))*sin(theta/2d0)
       sr(2,1,iop)=(-ii*nvec(1)+nvec(2))*sin(theta/2d0)
       sr(2,2,iop)=cos(theta/2d0)+ii*nvec(3)*sin(theta/2d0)

    end do

    do iop=1,nsymq
       sr(:,:,iop+nsymq)=-sr(:,:,iop)
    end do

    deallocate(a); deallocate(rwork); deallocate(w)

    ! output spin rotation matrices
    OPEN(11,FILE=trim(data_dir)//"/output/srg.dat",status="replace")

    do iop=1,nsymq
       write(11,*) iop
       write(11,'(4f10.6)') sr(1,:,iop)
       write(11,'(4f10.6)') sr(2,:,iop)
       write(11,*)""
    end do
    
    CLOSE(11)
    
    OPEN(12,FILE=trim(data_dir)//"/output/srg_import.txt",status="replace")

    write(12,fmt='(a)',advance='no')"{"
    do iop=1,nsymq
       write(12,'("{{",1f12.8,"+(",1f12.8," I),",1f12.8,"+(",1f12.8," I)},&
            & {",1f12.8,"+(",1f12.8," I),",1f12.8,"+(",1f12.8," I)}}")'&
            & ,advance='no') sr(1,:,iop),sr(2,:,iop)
       if(iop/=nsymq) write(12,fmt='(a)',advance='no')","
    end do
    write(12,fmt='(a)',advance='no')"}"
    
    CLOSE(12)

  end subroutine spin_rot

  ! specify the number of bands for calculation
  ! here, all bands in input files are cosidered
  subroutine get_occ_bands_full(NBAND,nband_occ)
    implicit none
    integer,intent(in) :: NBAND ! number of bands in input files
    integer,intent(out) :: nband_occ(Nk_irr) ! number of considered bands
    
    nband_occ=NBAND
    print*,"================="
    print*,""
    print*,"NBAND:",NBAND

  end subroutine get_occ_bands_full

  ! get how degeneracy occurs among considered bands 
  subroutine get_degeneracy_full_levels(Nk_irr,E_EIGI,NBAND,degeneracy,energy_level)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-point
    integer,intent(in) :: NBAND ! number of bands in input files
    real(8),intent(in) :: E_EIGI(NBAND,Nk_irr) ! energy of band
    integer,intent(out) :: degeneracy(NBAND,Nk_irr) ! number of degenerated bands
                                                    ! for each energy level
    integer,intent(out) :: energy_level(Nk_irr) ! number of energy level
                                                ! for each k-point

    integer :: ik,ib,count
    real(8),parameter :: err=1d-6

    degeneracy=0
    
    ! get degeneracy and number or energy levels
    do ik=1,Nk_irr
       count=1
       do ib=1,NBAND-1
             degeneracy(count,ik) = degeneracy(count,ik)+1
          if(E_EIGI(ib+1,ik) > E_EIGI(ib,ik)+err)then
             count = count+1
          end if
       end do
       degeneracy(count,ik) = degeneracy(count,ik)+1
       energy_level(ik) = count
    end do

  end subroutine get_degeneracy_full_levels

  ! get k-group G_k and the phase factor of represenatations  
  subroutine get_kgroup(Nk_irr,nsymq,SKI,rg,pg,nnp,gk,Gshift,kphase)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operation
    integer,intent(in) :: rg(3,3,nsymq) ! point group part of symmetry operation
                                        ! in the primitive coordinate
                                        ! of the reciprocal space
    integer,intent(in) :: pg(3,nsymq) ! translation part of symmetry operation
    integer,intent(in) :: nnp ! scale of fractional translation
    real(8),intent(in) :: SKI(3,Nk_irr) ! k-points
    logical,intent(out) :: gk(Nk_irr,nsymq) ! whether the symmetry operator
                                            ! is included in k-group
    integer,intent(out) :: Gshift(3,Nk_irr,nsymq) ! shift of k-vector caused
                                                  ! by the symmetry operation
    complex(8),intent(out) :: kphase(Nk_irr,nsymq) ! the phase factor of representations
                                                   ! in form of Exp(-ik t_g)
    
    integer :: ik,iop
    logical :: flag
    real(8) :: err=0.00001, tmp(3)
    real(8),parameter :: pi = acos(-1d0)
    complex(8),parameter :: ii = (0,1)
    
    gk = .false.
    flag = .false.
    kphase=1d0


    do ik=1,Nk_irr
       do iop=1,nsymq
          ! calculate the shift of k-vector caused by the symmetry operation
          tmp(1)=sum(rg(1,:,iop)*SKI(:,ik))-SKI(1,ik)
          tmp(2)=sum(rg(2,:,iop)*SKI(:,ik))-SKI(2,ik)
          tmp(3)=sum(rg(3,:,iop)*SKI(:,ik))-SKI(3,ik)
          Gshift(1,ik,iop)=nint(tmp(1))
          Gshift(2,ik,iop)=nint(tmp(2))
          Gshift(3,ik,iop)=nint(tmp(3))

          ! if the shift of k-vector is reciprocal lattice vector
          ! (integer in the basis of reciprocal lattice vector),
          ! that symmetry operation should be the element of k-group
          gk(ik,iop)=(abs(tmp(1)-Gshift(1,ik,iop))<err).and.&
               &(abs(tmp(2)-Gshift(2,ik,iop))<err).and.&
               &(abs(tmp(3)-Gshift(3,ik,iop))<err)

       end do       
    end do

    ! calculate kphase which is the translational part of representation
    do ik=1,Nk_irr
       do iop=1,nsymq
          kphase(ik,iop)=exp(-ii*sum(SKI(:,ik)*pg(:,iop))*2*pi/nnp)
       end do
    end do
       
  end subroutine get_kgroup

    
  ! get map of expansion reciprocal lattice vectors by symmetry operation
  ! and phase factor of representations
  subroutine get_transformed_G(Nk_irr,nsymq,maxNGp,NG_for_psi,KGI,gk,Gshift,pg,nnp,tfkg,phase)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: maxNGp ! maximum number among NG_for_psi
    integer,intent(in) :: NG_for_psi(Nk_irr) ! number of expansion reciprocal lattice vectors
    integer,intent(in) :: KGI(3,maxNGp,Nk_irr) ! expansion reciprocal lattice vectors
    integer,intent(in) :: Gshift(3,Nk_irr,nsymq) ! shift of k-vector caused
                                                 ! by the symmetry operation
    integer,intent(in) :: pg(3,nsymq) ! translation part of symmetry operation
    integer,intent(in) :: nnp ! scale of fractional translation
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether symmetry operation is inculded in k-group
    integer,intent(out) :: tfkg(maxNGp,nk_irr,nsymq) ! indices of expansion reciprocal lattice 
                                                     ! vectors transformed by symmetry operation 
    complex(8),intent(out) :: phase(maxNGp,Nk_irr,nsymq) ! the phase factor of representation
                                                         ! in form of Exp(-iG t_g) 

    integer :: ik,iop,ig
    real(8),parameter :: pi = acos(-1d0)
    complex(8),parameter :: ii = (0,1)

    tfkg=-1
    
    do ik=1,Nk_irr
       do iop=1,nsymq
          ! if considered operator is the element of k-group
          if(gk(ik,iop))then
             do ig=1,NG_for_psi(ik)

                ! calculate the phase factor of representation
                phase(ig,ik,iop)=exp(-ii*sum(KGI(:,ig,ik)*pg(:,iop))*2*pi/nnp)

                ! reciprocal lattice vector G' transformed by the symmetry operation
                tmpkg(1)=sum(rg(1,:,iop)*KGI(:,ig,ik))
                tmpkg(2)=sum(rg(2,:,iop)*KGI(:,ig,ik))
                tmpkg(3)=sum(rg(3,:,iop)*KGI(:,ig,ik))

                ! consider the shift of high-symmetry k-point 
                tmpkg(:)=tmpkg(:)+Gshift(:,ik,iop)

                ! search for index of transformed reciprocal lattice vector
                flag=.false.
                do i=1,NG_for_psi(ik)
                   if(tmpkg(1)==KGI(1,i,ik).and.&
                        tmpkg(2)==KGI(2,i,ik).and.&
                        tmpkg(3)==KGI(3,i,ik))then
                      tfkg(i,ik,iop)=ig
                      flag=.true.
                      exit
                   end if
                end do
                
             end do
          end if
       end do
    end do
  end subroutine get_transformed_G

  ! get and output spin factor systems
  subroutine get_spin_factor_system(Nk_irr,nsymq,tr,gk,sr,facsys_spin,product_list)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    real(8),intent(in) :: tr(3,3,nsymq) ! point group part of symmetry operation
                                        ! in the conventional coordinate
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether the symmetry operator
                                           ! is included in k-group
    complex(8),intent(in) :: sr(2,2,2*nsymq) ! spin rotation part of symmetry operators
    complex(8),intent(out) :: facsys_spin(Nk_irr,nsymq,nsymq) ! spin factor system
    integer,intent(out) :: product_list(nsymq,nsymq) ! product table of group operators
    
    integer :: ik,iop,jop,kop,i,j,mulop
    real(8) :: tmp(3,3),err=1d-5

    tmp=0d0
    product_list=0
    facsys_spin=0d0
    mulop=0

    
    do ik=1,Nk_irr
       do iop=1,nsymq
          if(gk(ik,iop))then
             do jop=1,nsymq
                ! if the operator is the element of k-group
                if(gk(ik,jop))then

                   tmp=0d0

                   ! search for the index k of the k-group element
                   ! g_k = g_i g_j
                   do i=1,3
                      do j=1,3
                         tmp(i,j)=sum(tr(i,:,iop)*tr(:,j,jop))
                      end do
                   end do

                   do kop=1,nsymq
                      if(maxval(abs(tmp(:,:)-tr(:,:,kop)))<err)then
                         mulop=kop
                         product_list(iop,jop)=kop
                         exit
                      end if
                   end do

                   ! calculate the spin factor system
                   if(abs(sr(1,1,mulop))>err)then
                      facsys_spin(ik,iop,jop)=sum(sr(1,:,iop)*sr(:,1,jop))/sr(1,1,mulop)
                   elseif(abs(sr(1,2,mulop))>err)then
                      facsys_spin(ik,iop,jop)=sum(sr(1,:,iop)*sr(:,2,jop))/sr(1,2,mulop)
                   else
                   end if
                end if
             end do
          end if
       end do
    end do
  end subroutine get_spin_factor_system
  
  ! get diagonal part of representation matrix
  subroutine representation_matrix(Nk_irr,nsymq,ncomp,maxNGI,NBAND,nband_occ,gk,&
       NGI,kphase,tfkg,phase,sr,C0,repmat)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: ncomp ! number of compornent of wavefunction
    integer,intent(in) :: maxNGI ! maximum number of NGI
    integer,intent(in) :: NBAND ! number of bands in input files
    integer,intent(in) :: nband_occ(Nk_irr) ! number of bands considered in calculation
    integer,intent(in) :: NGI(Nk_irr) ! expansion reciprocal lattice vector
    integer,intent(in) :: tfkg(maxNGI,Nk_irr,nsymq) ! expansion reciprocal lattice vector
                                                    ! transformed by symmetry operations
    complex(8),intent(in) :: sr(2,2,nsymq) ! spin rotation part of symmetry operation
    complex(8),intent(in) :: C0(ncomp,maxNGI,NBAND,Nk_irr) ! wavefunction
    complex(8),intent(in) :: kphase(Nk_irr,nsymq) ! phase factor in representation
    complex(8),intent(in) :: phase(maxNGI,Nk_irr,nsymq) ! phase factor in representation
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether the symmetry operator
                                           ! is included in k-group

    complex(8),intent(out) :: repmat(NBAND,Nk_irr,nsymq) ! diagonal part of representation matrix

    complex(8) :: tmp(2,maxNGI)

    integer :: ik, iop, ib

    repmat=0


    do ik=1,Nk_irr       
       do iop=1,nsymq
          if(gk(ik,iop))then
             do ib=1,nband_occ(ik)
                tmp=0d0
                
                ! if spin-orbit coupling is not considered in calculations
                if(ncomp==1)then
                   repmat(ib,ik,iop)=sum(conjg(C0(1,1:NGI(ik),ib,ik))*C0(1,tfkg(1:NGI(ik),ik,iop),ib,ik)&
                        *phase(1:NGI(ik),ik,iop))*kphase(ik,iop)
                   
                   ! if spin-orbit coupling is considered in calculations
                else if(ncomp==2)then
                   tmp(1,1:NGI(ik))=(sr(1,1,iop)*C0(1,tfkg(1:NGI(ik),ik,iop),ib,ik) &
                        + sr(1,2,iop)*C0(2,tfkg(1:NGI(ik),ik,iop),ib,ik)) &
                        *phase(1:NGI(ik),ik,iop)*kphase(ik,iop)
                   tmp(2,1:NGI(ik))=(sr(2,1,iop)*C0(1,tfkg(1:NGI(ik),ik,iop),ib,ik) &
                        + sr(2,2,iop)*C0(2,tfkg(1:NGI(ik),ik,iop),ib,ik)) &
                        *phase(1:NGI(ik),ik,iop)*kphase(ik,iop)
                   
                   repmat(ib,ik,iop)=sum(conjg(C0(1,1:NGI(ik),ib,ik))*tmp(1,1:NGI(ik)))+&
                        sum(conjg(C0(2,1:NGI(ik),ib,ik))*tmp(2,1:NGI(ik)))
                end if
                
             end do
          end if
       end do
    end do
    
    
  end subroutine representation_matrix

  ! output spin factor system  
  subroutine factor_system_out(Nk_irr,nsymq,gk,SKI,facsys_spin,data_dir)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether the symmetry operation
                                           ! is included in k-group
    real(8),intent(in) :: SKI(3,Nk_irr) ! k-points
    complex(8),intent(in) :: facsys_spin(Nk_irr,nsymq,nsymq) ! spin factor system
    character(150),intent(in) :: data_dir ! output directory
    
    integer :: ik,iop,jop
    
    ! output spin factor system
    OPEN(23,FILE=trim(data_dir)//"/output/factor_system_spin_import.txt",status="replace")
    
    write(23,fmt='(a)',advance='no')"{"

    do ik=1,Nk_irr
       write(23,fmt='("{{",1f8.5,",",1f8.5,",",1f8.5,"},")',advance='no')SKI(:,ik)
       write(23,fmt='(a)',advance='no')"{"
       do iop=1,nsymq
          write(23,fmt='(a)',advance='no')"{"
          do jop=1,nsymq
             if(gk(ik,iop).and.gk(ik,jop))then
                write(23,fmt='(1f8.5,"+(",1f8.5,")I")',advance='no')facsys_spin(ik,iop,jop)
             else
                write(23,fmt='("Null")',advance='no')
             end if
             if(jop/=nsymq)then
                write(23,fmt='(a)',advance='no')","
             end if
          end do
          write(23,fmt='(a)',advance='no')"}"
          if(iop/=nsymq)then
             write(23,fmt='(a)',advance='no')","
          end if
       end do
       
       write(23,fmt='(a)',advance='no')"}}"
       if(ik/=Nk_irr)then
          write(23,fmt='(a)',advance='no')","
       end if
    end do
    write(23,fmt='(a)',advance='no')"}"
    
    close(23)

    OPEN(21,FILE=trim(data_dir)//"/output/factor_system_spin.dat",status="replace")

    do ik=1,Nk_irr
       write(21,fmt='("k=(",1f8.5,",",1f8.5,",",1f8.5,"),")')SKI(:,ik)
       write(21,*)
       do iop=1,nsymq
          do jop=1,nsymq
             write(21,fmt='(" ",1f8.5,"+(",1f8.5,")i ")',advance='no')facsys_spin(ik,iop,jop)
          end do
          write(21,*)
       end do
       write(21,*)
    end do

    close(21)

    
  end subroutine factor_system_out
  
  ! get character 
  subroutine get_character(Nk_irr,nsymq,NBAND,energy_level,degeneracy,repmat,chrct)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: NBAND ! number of bands in input files
    integer,intent(in) :: energy_level(Nk_irr) ! number of energy level
    integer,intent(in) :: degeneracy(NBAND,Nk_irr) ! number of degeneracy
    complex(8),intent(in) :: repmat(NBAND,Nk_irr,nsymq) ! diagonal part of representation matrix
    complex(8),intent(out) :: chrct(Nk_irr,nsymq,maxval(energy_level)) ! character

    integer :: ik, iop, j, k, count 

    chrct=0d0
    count = 0

    ! take summation of diagonal part of representation matrix
    ! for each degenerated band of the same energy level
    do ik=1,Nk_irr
       do iop=1,nsymq
          count = 0
          do j=1,energy_level(ik)
             do k=1,degeneracy(j,ik)
                chrct(ik,iop,j)=chrct(ik,iop,j)&
                     &+repmat(count+k,ik,iop)
             end do
             count = count + degeneracy(j,ik)
          end do
       end do
    end do

  end subroutine get_character

  ! output character
  subroutine character_out(Nk_irr,nsymq,energy_level,SKI,gk,chrct,data_dir)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: energy_level(Nk_irr) ! number of energy level
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether the symmetry operator
                                           ! is included in k-group
    real(8),intent(in) :: SKI(3,Nk_irr) ! k-point
    complex(8),intent(in) :: chrct(Nk_irr,nsymq,maxval(energy_level)) ! character 
    character(150),intent(in) :: data_dir ! output directory
    
    integer :: ik,iop,ib

    
    ! output character
    OPEN(20,FILE=trim(data_dir)//"/output/character_import.txt",status="replace")

    
    write(20,fmt='(a)',advance='no')"{"

    do ik=1,Nk_irr
       write(20,fmt='("{{",1f8.5,",",1f8.5,",",1f8.5,"},")',advance='no')SKI(:,ik)
       write(20,fmt='(a)',advance='no')"{"
       do ib=1,energy_level(ik)
          write(20,fmt='(a)',advance='no')"{"

          do iop=1,nsymq
             if((gk(ik,iop).and.iop.le.nsymq))then
                write(20,fmt='(1f8.5,"+(",1f8.5,")I")',advance='no')chrct(ik,iop,ib)
             else
                write(20,fmt='("Null")',advance='no')
             end if
             if(iop/=nsymq)then
                write(20,fmt='(a)',advance='no')","
             end if
          end do

          write(20,fmt='(a)',advance='no')"}"
          if(ib/=energy_level(ik))then
             write(20,fmt='(a)',advance='no')","
          end if
       end do
       write(20,fmt='(a)',advance='no')"}}"
       if(ik/=Nk_irr)then
          write(20,fmt='(a)',advance='no')","
       end if
    end do
    write(20,fmt='(a)',advance='no')"}"
    
    CLOSE(20)

    OPEN(21,FILE=trim(data_dir)//"/output/character.dat",status="replace")

    do ik=1,Nk_irr
       write(21,fmt='("k=(",1f8.5,",",1f8.5,",",1f8.5,"),")')SKI(:,ik)
       write(21,*)"symmetry operator indices in G_k:"
       do iop=1,nsymq
          if(gk(ik,iop))then
             write(21,fmt='(" ",1I3," ")',advance='no')iop
          end if
       end do
       write(21,*)
       do ib=1,energy_level(ik)
          write(21,fmt='("level=",1I3," ")',advance='no')ib
          do iop=1,nsymq
             if(gk(ik,iop))then
                write(21,fmt='(" ",1f8.5,"+(",1f8.5,")i ")',advance='no')chrct(ik,iop,ib)
             else
                write(21,fmt='(" * ")',advance='no')
             end if
          end do
          write(21,*)
       end do
       write(21,*)
    end do
    
    CLOSE(21)

    
  end subroutine character_out
  
  ! get and output z4 index
  ! this process runs only when filling is specified
  subroutine get_z4(Nk_irr,nsymq,NBAND,iinv,filling,gk,energy_level,degeneracy,chrct,data_dir)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: NBAND ! number of bands in input files
    integer,intent(in) :: iinv ! index of inversion operator 
    integer,intent(in) :: filling ! number of filling
    integer,intent(in) :: energy_level(Nk_irr) ! number of energy level
    integer,intent(in) :: degeneracy(NBAND,Nk_irr) ! number of degeneracy
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether the symmetry operation
                                           ! is included in k-group
    complex(8),intent(in) :: chrct(Nk_irr,nsymq,maxval(energy_level)) ! character
    character(150),intent(in) :: data_dir ! output directory

    integer :: ik, il, count, occ, deno=4
    real(8) :: total, err
    integer :: z4
    logical :: err_flag

    total = 0d0
    count = 0
    occ = 0
    err=1d-6
    err_flag=.false.

    
    do ik=1,Nk_irr
       if(gk(ik,iinv))then
          count=count+1
       end if
    end do

    print*,"================="

    ! if inversion symmetry not exist
    if(iinv.lt.1)then
       print*,""
       print*,"Error:"
       print*,"There is no inversion symmetry. z4 index cannot be computed."
       err_flag = .true.

    ! if the number of specified TRIMs are less than 8
    else if(count.lt.8)then
 
       print*,""
       print*,"Error:"
       print*,"Only",count,"/ 8 TRIMs are found. Please specify 8 TRIMs."
       print*,"z4 index should be incorrect."
       err_flag = .true.

    ! if the number of specified TRIMs are greater than 8
    else if(count.gt.8)then

       print*,""
       print*,"Error:"
       print*,count,"/ 8 TRIMs are found. Please specify 8 different TRIMs."
       print*,"z4 index should be incorrect."
       err_flag = .true.
       
    end if

    ! sum of inversion parities of occupied bands
    do ik=1,Nk_irr
       if(gk(ik,iinv))then
          do il=1,energy_level(ik)
             if(filling<=sum(degeneracy(1:il,ik)))then
                occ=il
                exit
             end if
          end do
          
          total = total + sum(chrct(ik,iinv,1:occ))
       end if
    end do

    ! if the sum is not an integer number
    if(abs(total/deno-nint(total/deno))>err)then

       print*,""
       print*,"Error:"
       print*,"sum of parities is not an integer number"
       print*,"sum of parities:",total
       print*,total/deno
       print*,nint(total/deno)
       print*,abs(total/deno-nint(total/deno))
       print*,err
       err_flag = .true.
       
    end if

   
    z4=mod(nint(total/deno),4)

    if(z4<0) z4=z4+4

    ! if there is no error
    if(.not.err_flag) then
    
       ! output z4 index
       print*,""
       print*,"z4 index:",z4
    
       
       OPEN(20,FILE=trim(data_dir)//"/output/z4.dat",status="replace")
       
       write(20,*)"sum of parities for",count,"k-points:"
       write(20,*)total
       
       write(20,*)"z4 index:"
       write(20,*)z4
       
       CLOSE(20)

    end if

  end subroutine get_z4

  ! get and output character tables of irreducible representations for each k
  subroutine character_table_out(Nk_irr,nsymq,ncomp,energy_level,&
       rotA,nnp,pg,SKI,gk,chrct,product_list,facsys_spin,data_dir)
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: ncomp ! number of components of wavefunctions
    real(8),intent(in) :: rotA(3,3,nsymq) ! rotation part of symmetry operators
                                           ! with the primitive coorditnate
                                           ! of the real space
    integer,intent(in) :: pg(3,nsymq) ! translation part of symmetry operations
    integer,intent(in) :: nnp ! scale of fractional translations
    integer,intent(in) :: product_list(nsymq,nsymq) ! product table of group operators
    integer,intent(in) :: energy_level(Nk_irr) ! number of energy level
    logical,intent(inout) :: gk(Nk_irr,nsymq) ! whether the symmetry operation
                                           ! is included in k-group
    real(8),intent(in) :: SKI(3,Nk_irr) ! k-points
    complex(8),intent(in) :: chrct(Nk_irr,nsymq,maxval(energy_level)) ! character
    complex(8),intent(in) :: facsys_spin(Nk_irr,nsymq,nsymq) ! spin factor system
    character(150),intent(in) :: data_dir ! output directory

    real(8),parameter :: pi = acos(-1d0)
    complex(8),parameter :: ii = (0,1)
    real(8),parameter :: err=1d-6
    
    integer :: iv,ik,iop,jop,lop,ideg,c1,c2,c3,max_deg,irr_number,irr_level
    integer :: info, lda, lwork, n

    logical :: is_unique_chrct
    
    real(8) :: norm,vec(3)
    complex(8) :: n_irrep

    
    complex(8),allocatable :: a(:,:)
    real(8),allocatable ::  work(:)
    integer,allocatable :: ipiv(:)
    real(8),allocatable :: w(:)
    
    integer,allocatable :: irr_deg(:)
    complex(8),allocatable :: facsys_nonsym(:,:)
    complex(8),allocatable :: regrep(:,:,:)
    complex(8),allocatable :: regrep_inv(:,:,:)    
    real(8),allocatable :: seed_r(:,:)
    real(8),allocatable :: seed_i(:,:)
    complex(8),allocatable :: h_seed(:,:)
    complex(8),allocatable :: h_sym(:,:)
    complex(8),allocatable :: irr_diag(:,:)
    complex(8),allocatable :: irr_chrct(:,:)
    complex(8),allocatable :: irr_chrct_red(:,:)
    complex(8),allocatable :: tmp(:,:)
    
    OPEN(20,FILE=trim(data_dir)//"/output/irreps_list.dat",status="replace")
    OPEN(21,FILE=trim(data_dir)//"/output/irreps_number.dat",status="replace")
    
    
    do ik =1,Nk_irr

       ! order of G_k
       n=count(gk(ik,:))
       
       lda=n
       lwork=64*n
             
       allocate(irr_deg(n))
       allocate(facsys_nonsym(n,n))
       allocate(regrep(n,n,n))
       allocate(regrep_inv(n,n,n))
       allocate(seed_r(n,n))
       allocate(seed_i(n,n))
       allocate(h_seed(n,n))
       allocate(h_sym(n,n))
       allocate(tmp(n,n))

       allocate(irr_diag(n,n))
       
       allocate(a(n,n))
       allocate(ipiv(n))
       allocate(work(lwork))
       allocate(w(n))

       facsys_nonsym=0d0
       regrep=0d0
       h_seed=0d0
       h_sym=0d0
       
       c1=0
       c2=0

       ! get factor systems originate from fractional translation
       ! in nonsymmorphic space groups
       do iop = 1,nsymq
          if(gk(ik,iop))then
             c1=c1+1
             c2=0
             do jop = 1,nsymq
                if(gk(ik,jop))then
                   c2=c2+1
                   do iv=1,3
                      vec(iv)=pg(iv,jop)-sum(rotA(iv,:,iop)*pg(:,jop))
                   end do
                   facsys_nonsym(c1,c2)=exp(ii*2*pi*real(sum(SKI(:,ik)*vec(:)))/nnp)

                end if
             end do
          end if
       end do

       c1=0
       c2=0
       c3=0

       ! get regular representations
       do iop = 1,nsymq
          if(gk(ik,iop)) then
             c1=c1+1
             c2=0
             do jop = 1,nsymq
                if(gk(ik,jop)) then
                   c2=c2+1
                   c3=0
                   do lop = 1,nsymq
                      if(gk(ik,lop))then
                         c3=c3+1
                         if(product_list(lop,jop)==iop) then
                            if(ncomp==1)then
                               regrep(c1,c2,c3)=facsys_nonsym(c3,c2)
                            else
                               regrep(c1,c2,c3)=facsys_spin(ik,lop,jop)*facsys_nonsym(c3,c2)
                            end if
                         else
                            regrep(c1,c2,c3)=0d0
                         end if
                      end if
                   end do
                end if
             end do
          end if
       end do
       
       ! get inverse matrices of regular representations
       do lop=1,n

          a(:,:)=regrep(:,:,lop)

          call zgetrf(n,n,a,lda,ipiv,info)

          if(info==0)then
             call zgetri(n,a,lda,ipiv,work,lwork,info)
             regrep_inv(:,:,lop)=a(:,:)
          else
             print*,"Error for a regular representation matrix"
          end if
          
       end do
       
       ! get a seed of random hamiltonian
       call random_number(seed_r)
       call random_number(seed_i)

       ! get a random hermitian matrix
       do iop=1,n
          do jop=1,n
             h_seed(iop,jop)=seed_r(iop,jop)+seed_r(jop,iop)+ii*seed_i(iop,jop)-ii*seed_i(jop,iop)
          end do
       end do

       ! symmetrize the hermitian
       ! get a random hamiltonian with space group symmetry       
       tmp=0d0
       
       do lop=1,n
          do iop=1,n
             do jop=1,n
                tmp(iop,jop)=sum(h_seed(iop,:)*regrep_inv(:,jop,lop))
             end do
          end do
          do iop=1,n
             do jop=1,n
                h_sym(iop,jop)=h_sym(iop,jop)+sum(regrep(iop,:,lop)*tmp(:,jop))
             end do
          end do
       end do
       
       ! diagonalize the random hamiltonian
       a(:,:)=h_sym(:,:)

       call diag(n,a,w) 
       
       ! get degeneracy and number or energy levels for irreducible representations       
       c1=1
       irr_deg=0
       do ib=1,n-1
          irr_deg(c1) = irr_deg(c1)+1
          if(w(ib+1) > w(ib)+err)then
             c1 = c1+1
          end if
       end do
       irr_deg(c1) = irr_deg(c1)+1
       irr_level = c1
       
       ! normalize eigenvectors
       do iop=1,n
          norm=sqrt(sum(conjg(a(:,iop))*a(:,iop)))
          a(:,iop)=a(:,iop)/norm
       end do

       ! get the diagonal part of the representation matrix
       c1=0
       do lop=1,nsymq
          if(gk(ik,lop))then
             c1=c1+1
             do iop=1,n
                do jop=1,n
                   tmp(iop,jop)=sum(regrep(iop,:,c1)*a(:,jop))*exp(-ii*2*pi*sum(SKI(:,ik)*pg(:,lop)/nnp))
                end do
             end do
             do jop=1,n
                irr_diag(jop,c1)=sum(conjg(a(:,jop))*tmp(:,jop))
             end do
          end if
       end do

       ! get the characters of G_k
       allocate(irr_chrct(irr_level,n))
       allocate(irr_chrct_red(irr_level,n))

       irr_chrct=0d0

       c1=0
       do iop=1,irr_level
          do lop=1,irr_deg(iop)
             do jop=1,n
                irr_chrct(iop,jop)=irr_chrct(iop,jop)+irr_diag(c1+lop,jop)
             end do
          end do
          c1 = c1 + irr_deg(iop)
       end do

       ! max value of the degeneracy
       max_deg=maxval(irr_deg(:))

       ! sort by degeneracy
       ! delete duplicates
       irr_chrct_red=0d0
       
       c1=0
       c2=0
       do ideg=1,max_deg
          c2=0
          do iop=1,irr_level
             if(abs(irr_chrct(iop,1)-1d0*ideg)<err)then
                if(ideg==1 .or. c2==0) then
                   c1=c1+1
                   c2=c2+1
                   irr_chrct_red(c1,:)=irr_chrct(iop,:)
                else
                   is_unique_chrct=.true.
                   do jop=1,c2
                      if(maxval(abs(irr_chrct_red(c1-jop+1,:)-irr_chrct(iop,:)))<err) then
                         is_unique_chrct=.false.
                         exit
                      end if
                   end do
                   if(is_unique_chrct)then
                      c1=c1+1
                      c2=c2+1
                      irr_chrct_red(c1,:)=irr_chrct(iop,:)
                   end if
                end if
             end if
          end do
       end do

       ! number of irreducible representations
       irr_number=c1

       ! output the character tables of irreducible representations
       write(20,fmt='("k=(",1f8.5,",",1f8.5,",",1f8.5,"),")')SKI(:,ik)
       write(20,*)"symmetry operator indices in G_k:"
       do iop=1,nsymq
          if(gk(ik,iop))then
             write(20,fmt='(" ",1I3," ")',advance='no')iop
          end if
       end do
       write(20,*)
       do iop=1,irr_number
          write(20,fmt='("irrep ",1I3)',advance='no')iop
          do jop=1,n
             write(20,fmt='(" ",1f8.5,"+(",1f8.5,")i ")',advance='no')irr_chrct_red(iop,jop)
          end do
          write(20,*)
       end do
       write(20,*)
       
       ! output the number of irreducible representations for each energy level
       write(21,fmt='("k=(",1f8.5,",",1f8.5,",",1f8.5,"),")')SKI(:,ik)
       write(21,fmt='("number of irreducible representations:",1I3)')irr_number
       
       do ib=1,energy_level(ik)
          
          write(21,fmt='("level=",1I3," ")',advance='no')ib
          
          do iop=1,irr_number
             n_irrep=0d0
             c1=0
             do jop=1,nsymq
                if(gk(ik,jop))then
                   c1=c1+1
                   n_irrep=n_irrep+conjg(irr_chrct_red(iop,c1))*chrct(ik,jop,ib)
                end if
             end do
             write(21,fmt='(" ",1I3)',advance='no')nint(real(n_irrep)/n)             
             if(abs(nint(real(n_irrep)/n)-n_irrep/n)>err)then
                print*,"Warning: number of irreps is not an integer"
                print*,"at k=",SKI(:,ik),"and level=",ib,":"
                print*,n_irrep/n
             end if
          end do
          write(21,*)
       end do
       write(21,*)
       
       deallocate(irr_deg)
       deallocate(facsys_nonsym)
       deallocate(regrep)
       deallocate(regrep_inv)
       deallocate(seed_r)
       deallocate(seed_i)
       deallocate(h_seed)
       deallocate(h_sym)
       deallocate(tmp)

       deallocate(irr_diag)
       deallocate(irr_chrct)
       deallocate(irr_chrct_red)
       
       deallocate(a)
       deallocate(ipiv)
       deallocate(work)
       deallocate(w)

    end do
    
    CLOSE(20); CLOSE(21)
    
  end subroutine character_table_out

  
  ! output the file for CheckTopologicalMat
  ! this process runs only when filling is specified  
  subroutine Bilbao_out(Nk_irr,nsymq,ncomp,NBAND,filling,energy_level,&
       degeneracy,E_EIGI,nnp,pg,sr,rotA,SKI,gk,chrct,data_dir)
    implicit none
    integer,intent(in) :: Nk_irr ! number of k-points
    integer,intent(in) :: nsymq ! number of symmetry operations
    integer,intent(in) :: ncomp ! number of components of wavefunctions
    integer,intent(in) :: NBAND ! number of bands in input files
    integer,intent(in) :: filling ! number of filling
    integer,intent(in) :: energy_level(Nk_irr) ! number of energy levels
    integer,intent(in) :: degeneracy(NBAND,Nk_irr) ! number of degeneracy
    real(8),intent(in) :: rotA(3,3,nsymq) ! point group part of symmetry operation
                                          ! in the primitive coordinate
                                          ! of the real space
    integer,intent(in) :: pg(3,nsymq) ! translation part of symmetry operations
    integer,intent(in) :: nnp ! scale of fractional translations
    logical,intent(in) :: gk(Nk_irr,nsymq) ! whether the symmetry operation
                                           ! is included in k-group
    real(8),intent(in) :: SKI(3,Nk_irr) ! k-points
    real(8),intent(in) :: E_EIGI(NBAND,Nk_irr) ! energy of bands
    complex(8),intent(in) :: sr(2,2,nsymq) ! spin rotation part of symmetry
    complex(8),intent(in) :: chrct(Nk_irr,nsymq,maxval(energy_level)) ! character
    character(150),intent(in) :: data_dir ! output directory

    integer iop,ik,count
    real(8),parameter :: au_ev=27.2114 ! the constant to change the energy unit
                                       ! from a.u. to eV
    
    print*,"run Bilbao_out"
    
    ! generate the input file for CheckTopologicalMat
    OPEN(21,FILE=trim(data_dir)//"/output/Bilbao.txt",status="replace")

    write(21,'(1I3)') filling
    
    write(21,'(1I3)') ncomp-1
    write(21,'(1I3)') nsymq

    do iop=1,nsymq

       write(21,fmt='(3I3,3I3,3I3)',advance='no')nint(rotA(1,:,iop)),nint(rotA(2,:,iop)),nint(rotA(3,:,iop))

       
       write(21,fmt='(3f10.6)',advance='no')dble(pg(:,iop))/dble(nnp)

       write(21,fmt='(2f10.6)',advance='no')real(sr(1,1,iop)),aimag(sr(1,1,iop))
       write(21,fmt='(2f10.6)',advance='no')real(sr(1,2,iop)),aimag(sr(1,2,iop))
       write(21,fmt='(2f10.6)',advance='no')real(sr(2,1,iop)),aimag(sr(2,1,iop))
       write(21,fmt='(2f10.6)',advance='no')real(sr(2,2,iop)),aimag(sr(2,2,iop))


       write(21,*)
    end do
    write(21,'(1I3)') Nk_irr
    do ik=1,Nk_irr
       write(21,fmt='(3f10.6)')SKI(:,ik)
    end do

    
    do ik=1,Nk_irr
       count = 0
       do iop=1,nsymq
          if(gk(ik,iop)) count=count+1
       end do
       write(21,'(1I3)')count
       
       do iop=1,nsymq
          if(gk(ik,iop))then
             write(21,fmt='(1I3)',advance='no')iop
          end if
       end do
       write(21,*)

       count=1
       do ib=1,energy_level(ik)
          if(count>filling) exit
          write(21,fmt='(1I3)',advance='no')count
          write(21,fmt='(1I3)',advance='no')degeneracy(ib,ik)

          ! changing the unit of energy from a.u. to eV
          write(21,fmt='(1f11.6)',advance='no')E_EIGI(count,ik)*au_ev 
          
          count=count+degeneracy(ib,ik)
          
          do iop=1,nsymq
             if(gk(ik,iop))then
                write(21,fmt='(2f10.6)',advance='no')chrct(ik,iop,ib)
             end if
          end do
          write(21,*)
          
       end do
    end do

    CLOSE(21)
    
  end subroutine Bilbao_out


  ! diagonalize matrix
  subroutine diag(n,mat,w) 
    implicit none
    integer , intent(in) :: n 
    complex(8) , intent(inout) :: mat(n,n)
    real(8) , intent(out) :: w(n)
    integer :: lwork 
    integer :: lrwork 
    integer :: liwork 
    real(8) , allocatable :: rwork(:)
    complex(8) , allocatable :: work(:)
    integer , allocatable :: iwork(:)
    integer :: info 
    
    lwork = 2*n + n**2 + 3
    lrwork = 1 + 5*n + 2*n**2 + 3
    liwork = 3 + 5*n + 3
    allocate( rwork(lrwork) )
    allocate( work(lwork) )
    allocate( iwork(liwork) )
    
    info = 0
    call zheevd('V','U',n,mat,n,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    if(info /= 0) then
       write(6,*) 'info (subroutine diag):' , info
       stop
    end if
    
    deallocate(rwork,work,iwork)
    return 
  end subroutine diag
  
end program Main
