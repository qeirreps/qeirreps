program Main
  implicit none
  integer :: cmdcount, Nk_irr, NBAND, ncomp, nsymq, nnp, ik, ib, jb, ic, ig, iop, i, j, maxNGp, maxNGI, tmpkg(3), iinv,filling
  real(8) :: Ecut_for_psi, FermiEnergy, Etot, a1(3), a2(3), a3(3)
  logical :: flag,flag_fill
  character(150) :: data_dir,data_fill
  integer, allocatable ::  NGI(:), KGI(:,:,:), rg(:,:,:), pg(:,:),NG_for_psi(:), tfkg(:,:,:), nband_occ(:), degeneracy(:,:), energy_level(:),Gshift(:,:,:)
  real(8), allocatable :: SKI(:,:), E_EIGI(:,:), tr(:,:,:)
  complex(8), allocatable :: CO(:,:,:,:),repmat(:,:,:,:),chrct(:,:,:), sr(:,:,:),phase(:,:,:),kphase(:,:),facsys(:,:,:),facsys_spin(:,:,:)
  logical, allocatable :: gk(:,:)

  maxNGp = 0
  maxNGI = 0


!============================================

  ! get data

!============================================

  cmdcount=command_argument_count()
  if(cmdcount == 0)then
     stop"please specify a directory"
  else if(cmdcount == 1)then
     flag_fill=.false.
     call get_command_argument(1,data_dir)
     print*,"not specify the filling"
     print*,"reading:"//trim(data_dir)
  else if(cmdcount == 2)then
     flag_fill=.true.
     call get_command_argument(1,data_dir)
     call get_command_argument(2,data_fill)
     read(data_fill,*)filling
     print*,"reading:"//trim(data_dir)
     print*,"filling:",filling
  else
     stop"please specify single directory and filling"
  end if
  

  OPEN(117,FILE=trim(data_dir)//"/dir-wfn/dat.bandcalc")
  read(117,*) Ecut_for_psi
  read(117,*) FermiEnergy
  read(117,*) Etot
  CLOSE(117)


  OPEN(101,FILE=trim(data_dir)//"/dir-wfn/dat.sample-k")
  read(101,*) Nk_irr 
  allocate(SKI(3,Nk_irr))
  do ik=1,Nk_irr
     read(101,*) (SKI(i,ik),i=1,3)
  enddo
  CLOSE(101)


  OPEN(111,FILE=trim(data_dir)//"/dir-wfn/dat.eigenvalue")
  rewind(111)
  read(111,*) NBAND
  allocate(E_EIGI(NBAND,Nk_irr))
  do ik=1,Nk_irr
     do ib=1,NBAND
        read(111,*) E_EIGI(ib,ik)
     end do
  end do
  CLOSE(111)


  OPEN(132,FILE=trim(data_dir)//"/dir-wfn/dat.nkm")
  rewind(132)
  allocate(NGI(Nk_irr))
  do ik=1,Nk_irr
     read(132,*) NGI(ik)
  end do
  CLOSE(132)


  OPEN(104,FILE=trim(data_dir)//"/dir-wfn/dat.kg")
  allocate(NG_for_psi(Nk_irr))
  do ik=1,Nk_irr
     read(104,*) NG_for_psi(ik)
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
        read(104,*) (KGI(i,ig,ik),i=1,3) 
     end do
  end do
  CLOSE(104)


  OPEN(102,FILE=trim(data_dir)//"/dir-wfn/dat.wfn",FORM='unformatted')
  rewind(102)
  read(102) ncomp 

  print*,"ncomp=",ncomp
  print*,""

  maxNGI = maxval(NGI)
  allocate(CO(ncomp,maxNGI,NBAND,Nk_irr))
  CO=0.0d0

  do ik=1,Nk_irr
     do ib=1,NBAND 
        do ic=1,ncomp
           read(102)(CO(ic,i,ib,ik),i=1,NGI(ik))
        end do
     end do
  end do
  CLOSE(102)


  OPEN(105,FILE=trim(data_dir)//"/dir-wfn/dat.lattice")
  read(105,*) a1(1), a1(2), a1(3)
  read(105,*) a2(1), a2(2), a2(3)
  read(105,*) a3(1), a3(2), a3(3)
  CLOSE(105)


  OPEN(100,FILE=trim(data_dir)//"/dir-wfn/dat.symmetry")
  read(100,*) nsymq
  read(100,*) nnp
  allocate(rg(3,3,nsymq))
  allocate(pg(3,nsymq))
  do iop=1,nsymq
     read(100,*) ((rg(i,j,iop),i=1,3),j=1,3)
     read(100,*) (pg(i,iop),i=1,3)
  end do
  CLOSE(100)

!============================================

  ! Main

!============================================


  call get_reci_vec(a1,a2,a3)
  
  call operation_type(nsymq,rg,pg,iinv)

  allocate(tr(3,3,nsymq))
  call rg_out(nsymq,rg,a1,a2,a3,tr,data_dir)
  call pg_out(nsymq,nnp,pg,data_dir)

  allocate(sr(2,2,2*nsymq))
  call spin_rot(nsymq,tr,sr,data_dir)

  allocate(nband_occ(Nk_irr))
!  call get_occ_bands(Nk_irr,NBAND,E_EIGI,FermiEnergy,nband_occ)
  call get_occ_bands_full(Nk_irr,NBAND,E_EIGI,FermiEnergy,nband_occ)

  allocate(energy_level(Nk_irr))
  allocate(degeneracy(NBAND,Nk_irr))
!  call get_degeneracy(Nk_irr,E_EIGI,NBAND,nband_occ,degeneracy,energy_level)
  call get_degeneracy_full_levels(Nk_irr,E_EIGI,NBAND,nband_occ,degeneracy,energy_level)


  allocate(gk(Nk_irr,nsymq))
  allocate(Gshift(3,Nk_irr,nsymq))
  allocate(kphase(Nk_irr,nsymq))
  call get_kgroup(Nk_irr,nsymq,SKI,rg,pg,nnp,gk,Gshift,kphase)


  allocate(tfkg(maxNGp,Nk_irr,nsymq))
  allocate(phase(maxNGp,Nk_irr,nsymq))
  call get_transformed_G(Nk_irr,nsymq,maxNGp,NG_for_psi,KGI,gk,Gshift,pg,nnp,tfkg,phase)

  allocate(repmat(NBAND,NBAND,Nk_irr,ncomp*nsymq))
  call representation_matrix(Nk_irr,nsymq,ncomp,maxNGI,NBAND,nband_occ,NGI,SKI,kphase,tfkg,phase,sr,CO,repmat)
!  call repmat_out(Nk_irr,nsymq,ncomp,maxNGI,NBAND,nband_occ,NGI,SKI,repmat)


  allocate(facsys(Nk_irr,nsymq,nsymq))
  allocate(facsys_spin(Nk_irr,nsymq,nsymq))
  call get_factor_system(Nk_irr,nsymq,SKI,tr,pg,nnp,gk,facsys)
  call get_spin_factor_system(Nk_irr,nsymq,tr,gk,sr,facsys_spin)
  call factor_system_out(Nk_irr,nsymq,gk,SKI,facsys,facsys_spin,data_dir)
  

  call kgroup_out(Nk_irr,nsymq,SKI,gk,data_dir)

  allocate(chrct(Nk_irr,ncomp*nsymq,maxval(energy_level)))
  call get_character(Nk_irr,nsymq,ncomp,NBAND,energy_level,degeneracy,repmat,chrct)
  call character_out(Nk_irr,nsymq,ncomp,NBAND,energy_level,degeneracy,SKI,gk,chrct,data_dir)

  if(flag_fill) call get_kappa1(Nk_irr,nsymq,ncomp,NBAND,iinv,filling,gk,energy_level,degeneracy,chrct)

  print*,"================="
  print*,""
  print*,"job done"

!============================================

  ! End

!============================================

  deallocate(NGI); deallocate(KGI); deallocate(rg); deallocate(pg)
  deallocate(SKI); deallocate(E_EIGI)
  deallocate(CO)
  deallocate(nband_occ)
  deallocate(gk); deallocate(Gshift); deallocate(kphase)
  deallocate(tr); deallocate(sr)
  deallocate(degeneracy); deallocate(energy_level)
  deallocate(facsys); deallocate(facsys_spin)
  deallocate(repmat); deallocate(chrct)


contains

!============================================

  ! Subroutines

!============================================


  subroutine get_reci_vec(a1,a2,a3)
    implicit none
    real(8),intent(in) :: a1(3),a2(3),a3(3)

    real(8),parameter :: pi=acos(-1d0)
    real(8) :: vol,b(3,3)

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
  
  subroutine operation_type(nsymq,rg,pg,iinv)
    implicit none
    integer,intent(in) :: nsymq,  rg(3,3,nsymq), pg(3,nsymq)
    integer,intent(out) :: iinv
    integer :: i, j, k, iop, jop, kop, count1, count2, n=3
    logical :: flag(nsymq)
    real(8) :: err, tmp1(3,3),tmp2(3,3), id(3,3), inv(3,3), det

    err = 0.00001
    id = 0d0
    inv = 0d0
    iinv=0
    

    id(1,1)=1d0; id(2,2)=1d0; id(3,3)=1d0
    inv(1,1)=-1d0; inv(2,2)=-1d0; inv(3,3)=-1d0

    print*,"================="
    print*,"symmetry operation type:"

    do iop=1,nsymq

       det=rg(1,1,iop)*rg(2,2,iop)*rg(3,3,iop)&
            &+rg(1,2,iop)*rg(2,3,iop)*rg(3,1,iop)&
            &+rg(1,3,iop)*rg(2,1,iop)*rg(3,2,iop)&
            &-rg(1,1,iop)*rg(2,3,iop)*rg(3,2,iop)&
            &-rg(1,2,iop)*rg(2,1,iop)*rg(3,3,iop)&
            &-rg(1,3,iop)*rg(2,2,iop)*rg(3,1,iop)

       if(maxval(abs(id(:,:)-rg(:,:,iop)))<err)then
          if(maxval(abs(pg(:,iop)))>0)then
             print*,iop,"translation"
          else
             print*,iop,"E"
          end if
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

                   tmp1(:,:)=tmp2(:,:)
                   do i=1,3
                      do j=1,3
                         tmp2(i,j)=sum(tmp1(i,:)*rg(:,j,iop))
                      end do
                   end do
                   
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


  subroutine rg_out(nsymq,rg,a1,a2,a3,tr,data_dir)
    implicit none
    integer,intent(in) :: nsymq, rg(3,3,nsymq)
    real(8),intent(in) :: a1(3),a2(3),a3(3)
    real(8),intent(out) ::  tr(3,3,nsymq)
    character(150),intent(in) :: data_dir

    real(8) :: lat(3,3), latinv(3,3), tmp1(3,3), tmp2(3,3)
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

    print*,"================="
    print*,"space group symmetry:"
    print*,""
    print*,"rotation part:"
    
    call dgetrf(n,n,a,lda,ipiv,info)
    
    if(info==0)then
       call dgetri(n,a,lda,ipiv,work,lwork,info)
       latinv(:,:)=a(:,:)

    else
       print*,"The factor U is singular"
    end if
    
    do iop=1,nsymq

       a(:,:)=rg(:,:,iop)
       call dgetrf(n,n,a,lda,ipiv,info)
       
       if(info==0)then
          call dgetri(n,a,lda,ipiv,work,lwork,info)
          tmp1(:,:)=a(:,:)
       else
          print*,"The factor U is singular"
       end if

       do i=1,3
          do j= 1,3
             tmp2(i,j)=sum(tmp1(:,i)*latinv(:,j))
          end do
       end do

       do i=1,3
          do j= 1,3
             tr(i,j,iop)=sum(lat(i,:)*tmp2(:,j))
          end do
       end do
       
    end do
    

    do iop=1,nsymq
       print*, iop
       print'(3f6.2)',tr(1,:,iop)
       print'(3f6.2)',tr(2,:,iop)
       print'(3f6.2)',tr(3,:,iop)
       print*,""
    end do


    OPEN(12,FILE=trim(data_dir)//"/output/rg_import.txt",status="replace")

    write(12,fmt='(a)',advance='no')"{"
    do iop=1,nsymq
       write(12,'("{{",1f12.8,",",1f12.8,",",1f12.8,"},{",1f12.8,",",1f12.8,",",1f12.8,"},{",1f12.8,",",1f12.8,",",1f12.8,"}}")',advance='no') tr(1,:,iop),tr(2,:,iop),tr(3,:,iop)
       if(iop/=nsymq)then
          write(12,fmt='(a)',advance='no')","
       end if
    end do
    write(12,fmt='(a)',advance='no')"}"
    
    CLOSE(12)

    deallocate(a,work,ipiv)

  end subroutine rg_out


  subroutine pg_out(nsymq,nnp,pg,data_dir)
    implicit none
    integer,intent(in) :: nsymq, nnp, pg(3,nsymq)
    character(150),intent(in) :: data_dir
    
    integer :: i, j, iop

    print*,"translation part:"
    

    do iop=1,nsymq
       print*,iop
       print'(1I3,"/",1I3," ",1I3,"/",1I3," ",1I3,"/",1I3)',pg(1,iop),nnp,pg(2,iop),nnp,pg(3,iop),nnp
       print*,""
    end do


    OPEN(12,FILE=trim(data_dir)//"/output/tg_import.txt",status="replace")

    write(12,fmt='(a)',advance='no')"{"
    do iop=1,nsymq
       write(12,'("{",1I3,"/",1I3,",",1I3,"/",1I3,",",1I3,"/",1I3,"}")',advance='no') pg(1,iop),nnp,pg(2,iop),nnp,pg(3,iop),nnp
       if(iop/=nsymq)then
          write(12,fmt='(a)',advance='no')","
       end if
    end do
    write(12,fmt='(a)',advance='no')"}"
    
    CLOSE(12)

  end subroutine pg_out


  subroutine spin_rot(nsymq,tr,sr,data_dir)
    implicit none
    integer,intent(in) :: nsymq
    real(8),intent(in) ::  tr(3,3,nsymq)
    complex(8),intent(out) ::  sr(2,2,2*nsymq)
    character(150),intent(in) :: data_dir

    complex,parameter :: ii=(0,1)
    integer :: i, j, k, iop, jop, kop, count1, count2, n=3, lda, lwork, info
    logical :: flag(nsymq)
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

    id(1,1)=1d0; id(2,2)=1d0; id(3,3)=1d0
    inv(1,1)=-1d0; inv(2,2)=-1d0; inv(3,3)=-1d0

    print*,'spin rotation part:'
    
    do iop=1,nsymq
       tmp1=0d0; tmp2=0d0
       vec=0d0; nvec=0d0; norm=0d0; trace=0d0; theta=0

       det=tr(1,1,iop)*tr(2,2,iop)*tr(3,3,iop)&
            &+tr(1,2,iop)*tr(2,3,iop)*tr(3,1,iop)&
            &+tr(1,3,iop)*tr(2,1,iop)*tr(3,2,iop)&
            &-tr(1,1,iop)*tr(2,3,iop)*tr(3,2,iop)&
            &-tr(1,2,iop)*tr(2,1,iop)*tr(3,3,iop)&
            &-tr(1,3,iop)*tr(2,2,iop)*tr(3,1,iop)

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
       
       if(maxval(abs(tmp1(:,:)-id(:,:)))<err)then
          sr(:,:,iop)=reshape((/1d0,0d0,0d0,1d0/),(/2,2/))
          cycle
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
       else
          vec(1)=tmp2(3,2)
          vec(2)=tmp2(1,3)
          vec(3)=tmp2(2,1)
       end if

       
       norm=sqrt(sum(vec(:)**2))
       nvec(:)=vec(:)/norm
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
       
       sr(1,1,iop)=cos(theta/2d0)-ii*nvec(3)*sin(theta/2d0)
       sr(1,2,iop)=-(ii*nvec(1)+nvec(2))*sin(theta/2d0)
       sr(2,1,iop)=(-ii*nvec(1)+nvec(2))*sin(theta/2d0)
       sr(2,2,iop)=cos(theta/2d0)+ii*nvec(3)*sin(theta/2d0)

    end do

    do iop=1,nsymq
       sr(:,:,iop+nsymq)=-sr(:,:,iop)
    end do

    deallocate(a); deallocate(rwork); deallocate(w)


    do iop=1,nsymq
       print*, iop
       print'(4f10.6)', sr(1,:,iop)
       print'(4f10.6)', sr(2,:,iop)
       print*,""
    end do

    OPEN(12,FILE=trim(data_dir)//"/output/srg_import.txt",status="replace")

    write(12,fmt='(a)',advance='no')"{"
    do iop=1,nsymq
       write(12,'("{{",1f12.8,"+(",1f12.8," I),",1f12.8,"+(",1f12.8," I)},{",1f12.8,"+(",1f12.8," I),",1f12.8,"+(",1f12.8," I)}}")',advance='no') sr(1,:,iop),sr(2,:,iop)
       if(iop/=nsymq) write(12,fmt='(a)',advance='no')","
    end do
    write(12,fmt='(a)',advance='no')"}"
    
    CLOSE(12)

  end subroutine spin_rot

  
  subroutine get_occ_bands(Nk_irr,NBAND,E_EIGI,FermiEnergy,nband_occ)
    implicit none
    integer,intent(in) :: Nk_irr, NBAND
    real(8),intent(in) :: E_EIGI(NBAND,Nk_irr), FermiEnergy
    integer,intent(out) :: nband_occ(Nk_irr)
    
    integer :: ik,ib

    nband_occ=0
    
    do ik=1,Nk_irr
       do ib=1,NBAND
          if(E_EIGI(ib,ik) > FermiEnergy)then
             nband_occ(ik)=ib-1
             print*,FermiEnergy,E_EIGI(ib-1,ik),E_EIGI(ib,ik),ib-1
             exit
          end if
       end do
    end do
  end subroutine get_occ_bands


  subroutine get_occ_bands_full(Nk_irr,NBAND,E_EIGI,FermiEnergy,nband_occ)
    implicit none
    integer,intent(in) :: Nk_irr, NBAND
    real(8),intent(in) :: E_EIGI(NBAND,Nk_irr), FermiEnergy
    integer,intent(out) :: nband_occ(Nk_irr)
    
    integer :: ik,ib

    nband_occ=NBAND
    print*,"================="
    print*,""
    print*,"NBAND:",NBAND

  end subroutine get_occ_bands_full

 
  subroutine get_degeneracy(Nk_irr,E_EIGI,NBAND,nband_occ,degeneracy,energy_level)
    implicit none
    integer,intent(in) :: Nk_irr,NBAND,nband_occ(Nk_irr)
    real(8),intent(in) :: E_EIGI(NBAND,Nk_irr)
    integer,intent(out) :: degeneracy(NBAND,Nk_irr), energy_level(Nk_irr)

    integer :: ik,ib,count
    real(8),parameter :: err=1d-10

    degeneracy=0
    
    do ik=1,Nk_irr
       count=1
       do ib=1,nband_occ(ik)
             degeneracy(count,ik) = degeneracy(count,ik)+1
          if(E_EIGI(ib+1,ik) > E_EIGI(ib,ik)+err)then
             count = count+1
          end if
       end do
       energy_level(ik) = count-1
!       print*,ik,(sum(degeneracy(:,ik))==nband_occ(ik)),sum(degeneracy(:,ik)),nband_occ(ik)
    end do

  end subroutine get_degeneracy


  subroutine get_degeneracy_full_levels(Nk_irr,E_EIGI,NBAND,nband_occ,degeneracy,energy_level)
    implicit none
    integer,intent(in) :: Nk_irr,NBAND,nband_occ(Nk_irr)
    real(8),intent(in) :: E_EIGI(NBAND,Nk_irr)
    integer,intent(out) :: degeneracy(NBAND,Nk_irr), energy_level(Nk_irr)

    integer :: ik,ib,count
    real(8),parameter :: err=1d-6

    degeneracy=0
    
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

  
  subroutine get_kgroup(Nk_irr,nsymq,SKI,rg,pg,nnp,gk,Gshift,kphase)
    implicit none
    integer,intent(in) :: Nk_irr, nsymq,rg(3,3,nsymq),pg(3,nsymq),nnp
    real(8),intent(in) :: SKI(3,Nk_irr)
    logical,intent(out) :: gk(Nk_irr,nsymq)
    integer,intent(out) :: Gshift(3,Nk_irr,nsymq)
    complex(8),intent(out) :: kphase(Nk_irr,nsymq)
    
    integer :: ik,iop,jop,kop,tmp2(3)
    logical :: flag
    real(8) :: err=0.00001, tmp(3)
    real(8),parameter :: pi = acos(-1d0)
    complex(8),parameter :: ii = (0,1)
    
    gk = .false.
    flag = .false.
    kphase=1d0


    do ik=1,Nk_irr
       do iop=1,nsymq
          tmp(1)=sum(rg(1,:,iop)*SKI(:,ik))-SKI(1,ik)
          tmp(2)=sum(rg(2,:,iop)*SKI(:,ik))-SKI(2,ik)
          tmp(3)=sum(rg(3,:,iop)*SKI(:,ik))-SKI(3,ik)
          Gshift(1,ik,iop)=nint(tmp(1))
          Gshift(2,ik,iop)=nint(tmp(2))
          Gshift(3,ik,iop)=nint(tmp(3))

          gk(ik,iop)=(abs(tmp(1)-Gshift(1,ik,iop))<err).and.&
               &(abs(tmp(2)-Gshift(2,ik,iop))<err).and.&
               &(abs(tmp(3)-Gshift(3,ik,iop))<err)

       end do

       !option for calculating representation of G_k/T, not G_k
       !from here
       
!       flag=.false.

!       do iop=1,nsymq
!          do jop=1,nsymq

!             tmp2(1)=sum(rg(1,:,iop)*pg(:,jop))
!             tmp2(2)=sum(rg(2,:,iop)*pg(:,jop))
!             tmp2(3)=sum(rg(3,:,iop)*pg(:,jop))

!             if(gk(ik,iop).and.gk(ik,jop).and.mod(sum(Gshift(:,ik,iop)*tmp2(:)),nnp)/=0)then
!                print*,"k:",ik
!                do kop=1,nsymq
!                   kphase(ik,kop)=exp(-ii*sum(SKI(:,ik)*pg(:,kop))*2*pi/nnp)
!                end do
!                flag=.true.
!                exit
!             end if
!          end do
!          if(flag)exit
!       end do

       !to here
       
    end do

    do ik=1,Nk_irr
       do iop=1,nsymq
          kphase(ik,iop)=exp(-ii*sum(SKI(:,ik)*pg(:,iop))*2*pi/nnp)
       end do
    end do
       
  end subroutine get_kgroup

    
  subroutine get_transformed_G(Nk_irr,nsymq,maxNGp,NG_for_psi,KGI,gk,Gshift,pg,nnp,tfkg,phase)
    implicit none
    integer,intent(in) :: Nk_irr,nsymq,maxNGp,NG_for_psi(Nk_irr),KGI(3,maxNGp,Nk_irr),Gshift(3,Nk_irr,nsymq), pg(3,nsymq), nnp
    logical,intent(in) :: gk(Nk_irr,nsymq)
    integer,intent(out) :: tfkg(maxNGp,nk_irr,nsymq)
    complex(8),intent(out) :: phase(maxNGp,Nk_irr,nsymq)

    integer :: ik,iop,ig
    real(8),parameter :: pi = acos(-1d0)
    complex(8),parameter :: ii = (0,1)

    tfkg=-1
    
    do ik=1,Nk_irr
       do ig =1,NG_for_psi(ik)
       end do
    end do

    
    do ik=1,Nk_irr
       do iop=1,nsymq
          if(gk(ik,iop))then
             do ig=1,NG_for_psi(ik)
                
                phase(ig,ik,iop)=exp(-ii*sum(KGI(:,ig,ik)*pg(:,iop))*2*pi/nnp)
                
                tmpkg(1)=sum(rg(1,:,iop)*KGI(:,ig,ik))
                tmpkg(2)=sum(rg(2,:,iop)*KGI(:,ig,ik))
                tmpkg(3)=sum(rg(3,:,iop)*KGI(:,ig,ik))
                
                tmpkg(:)=tmpkg(:)+Gshift(:,ik,iop)
                
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
                if(.not.flag)then
                end if
                
             end do
          end if
       end do
    end do
  end subroutine get_transformed_G

  subroutine get_factor_system(Nk_irr,nsymq,SKI,tr,pg,nnp,gk,facsys)
    implicit none
    integer,intent(in) :: Nk_irr, nsymq,pg(3,nsymq),nnp
    real(8),intent(in) :: tr(3,3,nsymq),SKI(3,Nk_irr)
    logical,intent(in) :: gk(Nk_irr,nsymq)
    complex(8),intent(out) :: facsys(Nk_irr,nsymq,nsymq)
    
    integer :: ik,iop,jop
    real(8) :: tmp(3)
    real(8),parameter :: pi = acos(-1d0)
    complex(8),parameter :: ii = (0,1)
    
    facsys=0d0
    tmp=0d0
    
    do ik=1,Nk_irr
       do iop=1,nsymq
          if(gk(ik,iop))then
             do jop=1,nsymq
                if(gk(ik,jop))then
                   tmp(1)=sum(tr(1,:,iop)*pg(:,jop))
                   tmp(2)=sum(tr(2,:,iop)*pg(:,jop))
                   tmp(3)=sum(tr(3,:,iop)*pg(:,jop))
                   facsys(ik,iop,jop)=exp(-ii*sum(SKI(:,ik)*(tmp(:)-pg(:,jop)))*2*pi/nnp)
                end if
             end do
          end if
       end do
    end do
                   
  end subroutine get_factor_system

  subroutine get_spin_factor_system(Nk_irr,nsymq,tr,gk,sr,facsys_spin)
    implicit none
    integer,intent(in) :: Nk_irr, nsymq
    real(8),intent(in) :: tr(3,3,nsymq)
    logical,intent(in) :: gk(Nk_irr,nsymq)
    complex(8),intent(in) :: sr(2,2,2*nsymq)
    complex(8),intent(out) :: facsys_spin(Nk_irr,nsymq,nsymq)
    
    integer :: ik,iop,jop,kop,i,j,mulop
    real(8) :: tmp(3,3),err=1d-5

    tmp=0d0
    facsys_spin=0d0
    mulop=0

    
    do ik=1,Nk_irr
       do iop=1,nsymq
          if(gk(ik,iop))then
             do jop=1,nsymq
                if(gk(ik,jop))then

                   tmp=0d0

                   do i=1,3
                      do j=1,3
                         tmp(i,j)=sum(tr(i,:,iop)*tr(:,j,jop))
                      end do
                   end do

                   do kop=1,nsymq
                      if(maxval(abs(tmp(:,:)-tr(:,:,kop)))<err)then
                         mulop=kop
                         exit
                      end if
                   end do

                   if(abs(sr(1,1,mulop))>err)then
                      facsys_spin(ik,iop,jop)=sum(sr(1,:,iop)*sr(:,1,jop))/sr(1,1,mulop)
                   elseif(abs(sr(1,2,mulop))>err)then
                      facsys_spin(ik,iop,jop)=sum(sr(1,:,iop)*sr(:,2,jop))/sr(1,2,mulop)
                   else
!                      print*,"error on",iop,jop,mulop
                   end if
                end if
             end do
          end if
       end do
    end do
  end subroutine get_spin_factor_system

  
  subroutine representation_matrix(Nk_irr,nsymq,ncomp,maxNGI,NBAND,nband_occ,NGI,SKI,kphase,tfkg,phase,sr,CO,repmat)
    implicit none
    integer,intent(in) :: Nk_irr,nsymq,ncomp,maxNGI,NBAND,nband_occ(Nk_irr),NGI(Nk_irr),tfkg(maxNGI,Nk_irr,nsymq)
    real(8),intent(in) :: SKI(3,Nk_irr)
    complex(8),intent(in) ::  sr(2,2,nsymq),CO(ncomp,maxNGI,NBAND,Nk_irr),kphase(Nk_irr,nsymq),phase(maxNGI,Nk_irr,nsymq)
    complex(8),intent(out) :: repmat(NBAND,NBAND,Nk_irr,ncomp*nsymq)

    complex(8) :: tmp(2,maxNGI)

    integer :: ik, iop, ib, jb

    repmat=0

    
    do ik=1,Nk_irr       
       do iop=1,nsymq
          do ib=1,nband_occ(ik)
             do jb=1,nband_occ(ik)
                tmp=0d0

                if(ncomp==1)then
                   repmat(ib,jb,ik,iop)=sum(conjg(CO(1,1:NGI(ik),ib,ik))*CO(1,tfkg(1:NGI(ik),ik,iop),jb,ik)*phase(1:NGI(ik),ik,iop))*kphase(ik,iop)
                else if(ncomp==2)then
                   tmp(1,1:NGI(ik))=(sr(1,1,iop)*CO(1,tfkg(1:NGI(ik),ik,iop),jb,ik) + sr(1,2,iop)*CO(2,tfkg(1:NGI(ik),ik,iop),jb,ik))*phase(1:NGI(ik),ik,iop)*kphase(ik,iop)
                   tmp(2,1:NGI(ik))=(sr(2,1,iop)*CO(1,tfkg(1:NGI(ik),ik,iop),jb,ik) + sr(2,2,iop)*CO(2,tfkg(1:NGI(ik),ik,iop),jb,ik))*phase(1:NGI(ik),ik,iop)*kphase(ik,iop)

                   repmat(ib,jb,ik,iop)=sum(conjg(CO(1,1:NGI(ik),ib,ik))*tmp(1,1:NGI(ik)))+&
                        &sum(conjg(CO(2,1:NGI(ik),ib,ik))*tmp(2,1:NGI(ik)))
                   
                   tmp=0d0
                   tmp(1,1:NGI(ik))=(sr(1,1,iop+nsymq)*CO(1,tfkg(1:NGI(ik),ik,iop),jb,ik) + sr(1,2,iop+nsymq)*CO(2,tfkg(1:NGI(ik),ik,iop),jb,ik))*phase(1:NGI(ik),ik,iop)*kphase(ik,iop)
                   tmp(2,1:NGI(ik))=(sr(2,1,iop+nsymq)*CO(1,tfkg(1:NGI(ik),ik,iop),jb,ik) + sr(2,2,iop+nsymq)*CO(2,tfkg(1:NGI(ik),ik,iop),jb,ik))*phase(1:NGI(ik),ik,iop)*kphase(ik,iop)

                   repmat(ib,jb,ik,iop+nsymq)=sum(conjg(CO(1,1:NGI(ik),ib,ik))*tmp(1,1:NGI(ik)))+&
                        &sum(conjg(CO(2,1:NGI(ik),ib,ik))*tmp(2,1:NGI(ik)))
                end if

             end do
          end do
       end do
    end do
    
    
  end subroutine representation_matrix

  
  subroutine factor_system_out(Nk_irr,nsymq,gk,SKI,facsys,facsys_spin,data_dir)
    implicit none
    integer,intent(in) :: Nk_irr,nsymq
    logical,intent(in) :: gk(Nk_irr,nsymq)
    real(8),intent(in) :: SKI(3,Nk_irr)
    complex(8),intent(in) :: facsys(Nk_irr,nsymq,nsymq),facsys_spin(Nk_irr,nsymq,nsymq)
    character(150),intent(in) :: data_dir
    
    integer :: ik,iop,jop
    
    OPEN(20,FILE=trim(data_dir)//"/output/factor_system_nonsymmorphic_import.txt",status="replace")
    
    write(20,fmt='(a)',advance='no')"{"

    do ik=1,Nk_irr
       write(20,fmt='("{{",1f8.5,",",1f8.5,",",1f8.5,"},")',advance='no')SKI(:,ik)
       write(20,fmt='(a)',advance='no')"{"
       do iop=1,nsymq
          write(20,fmt='(a)',advance='no')"{"
          do jop=1,nsymq
             if(gk(ik,iop).and.gk(ik,jop))then
                write(20,fmt='(1f8.5,"+(",1f8.5,")I")',advance='no')facsys(ik,iop,jop)
             else
                write(20,fmt='("Null")',advance='no')
             end if
             if(jop/=nsymq)then
                write(20,fmt='(a)',advance='no')","
             end if
          end do
          write(20,fmt='(a)',advance='no')"}"
          if(iop/=nsymq)then
             write(20,fmt='(a)',advance='no')","
          end if
       end do
       
       write(20,fmt='(a)',advance='no')"}}"
       if(ik/=Nk_irr)then
          write(20,fmt='(a)',advance='no')","
       end if
    end do
    write(20,fmt='(a)',advance='no')"}"
    
    close(20)

    OPEN(21,FILE=trim(data_dir)//"/output/factor_system_nonsymmorphic.dat",status="replace")

    do ik=1,Nk_irr
       write(21,fmt='("k=(",1f8.5,",",1f8.5,",",1f8.5,"),")')SKI(:,ik)
       write(21,*)""
       do iop=1,nsymq
          do jop=1,nsymq
             write(21,fmt='(" ",1f8.5,"+(",1f8.5,")i ")',advance='no')facsys(ik,iop,jop)
          end do
          write(21,*)
       end do
       write(21,*)
    end do

    close(21)
    
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
  
  
  subroutine get_character(Nk_irr,nsymq,ncomp,NBAND,energy_level,degeneracy,repmat,chrct)
    implicit none
    integer,intent(in) :: Nk_irr, nsymq, ncomp, NBAND, energy_level(Nk_irr), degeneracy(NBAND,Nk_irr)
    complex(8),intent(in) :: repmat(NBAND,NBAND,Nk_irr,ncomp*nsymq)
    complex(8),intent(out) :: chrct(Nk_irr,ncomp*nsymq,maxval(energy_level))

    integer :: ik, iop, j, k, count 

    chrct=0d0
    count = 0

    do ik=1,Nk_irr
       do iop=1,ncomp*nsymq
          count = 0
          do j=1,energy_level(ik)
             do k=1,degeneracy(j,ik)
                chrct(ik,iop,j)=chrct(ik,iop,j)&
                     &+repmat(count+k,count+k,ik,iop)
             end do
             count = count + degeneracy(j,ik)
          end do
       end do
    end do

  end subroutine get_character

  
  subroutine get_kappa1(Nk_irr,nsymq,ncomp,NBAND,iinv,filling,gk,energy_level,degeneracy,chrct)
    implicit none
    integer,intent(in) :: Nk_irr, nsymq, ncomp, NBAND, iinv, filling, energy_level(Nk_irr), degeneracy(NBAND,Nk_irr)
    logical,intent(in) :: gk(Nk_irr,nsymq)
    complex(8),intent(in) :: chrct(Nk_irr,ncomp*nsymq,maxval(energy_level))

    integer :: ik, il, count, occ
    real(8) :: kappa1

    kappa1 = 0d0
    count = 0
    occ = 0

    do ik=1,Nk_irr
       if(gk(ik,iinv))then
          count = count+1
          do il=1,energy_level(ik)
             if(filling<=sum(degeneracy(1:il,ik)))then
                occ=il
                exit
             end if
          end do
          
          kappa1 = kappa1 + sum(chrct(ik,iinv,1:occ))
       end if
    end do

    print*,"================="
    print*,""
    print*,"summation of parities for",count,"k-points"
    print*,kappa1

  end subroutine get_kappa1

  

  subroutine character_out(Nk_irr,nsymq,ncomp,NBAND,energy_level,degeneracy,SKI,gk,chrct,data_dir)
    implicit none
    integer,intent(in) :: Nk_irr,nsymq,ncomp,NBAND,energy_level(Nk_irr),degeneracy(NBAND,Nk_irr)
    logical,intent(in) :: gk(Nk_irr,nsymq)
    real(8),intent(in) :: SKI(3,Nk_irr)
    complex(8),intent(in) :: chrct(Nk_irr,ncomp*nsymq,maxval(energy_level))
    character(150),intent(in) :: data_dir
    
    integer :: ik,iop,ib
    
    OPEN(20,FILE=trim(data_dir)//"/output/character_import.txt",status="replace")

    
    write(20,fmt='(a)',advance='no')"{"

    do ik=1,Nk_irr
       write(20,fmt='("{{",1f8.5,",",1f8.5,",",1f8.5,"},")',advance='no')SKI(:,ik)
       write(20,fmt='(a)',advance='no')"{"
       do ib=1,energy_level(ik)
          write(20,fmt='(a)',advance='no')"{"

          do iop=1,ncomp*nsymq
             if((gk(ik,iop).and.iop.le.nsymq).or.(gk(ik,iop-nsymq).and.iop.gt.nsymq))then
                write(20,fmt='(1f8.5,"+(",1f8.5,")I")',advance='no')chrct(ik,iop,ib)
             else
                write(20,fmt='("Null")',advance='no')
             end if
             if(iop/=ncomp*nsymq)then
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
    
    close(20)

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
                write(20,fmt='(" * ")',advance='no')
             end if
          end do
          write(21,*)
       end do
       write(21,*)
    end do
    
    CLOSE(21)

    
  end subroutine character_out


  subroutine kgroup_out(Nk_irr,nsymq,SKI,gk,data_dir)
    integer,intent(in) :: Nk_irr,nsymq
    logical,intent(in) :: gk(Nk_irr,nsymq)
    real(8),intent(in) :: SKI(3,Nk_irr)
    character(150),intent(in) :: data_dir

    integer :: ik,iop,lastsym

    OPEN(21,FILE=trim(data_dir)//"/output/kg_import.txt",status="replace")

    write(21,fmt='(a)',advance='no')"{"
    do ik=1,Nk_irr
       lastsym=1
       write(21,fmt='("{{",1f8.5,",",1f8.5,",",1f8.5,"},{")',advance='no')SKI(:,ik)
       do iop=1,nsymq
          if(gk(ik,iop))then
             lastsym=iop
          end if
       end do
       do iop=1,nsymq
          if(iop==lastsym)then
             write(21,fmt='(1I3,"}}")',advance='no')iop
             exit
          else if(gk(ik,iop))then
             write(21,fmt='(1I3,",")',advance='no')iop
          end if
       end do
       if(ik/=Nk_irr)then
          write(21,fmt='(a)',advance='no')","
       end if
    end do
    write(21,fmt='(a)',advance='no')"}"

    CLOSE(21)
    
  end subroutine kgroup_out


end program Main
