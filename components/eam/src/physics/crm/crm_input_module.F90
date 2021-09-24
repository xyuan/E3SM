module crm_input_module
   use params_kind, only: crm_rknd
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif
   use openacc_utils
   implicit none
   private
   public crm_input_type
   public crm_input_initialize
   public crm_input_finalize

#ifndef MODAL_AERO
   integer, parameter :: ntot_amode=1
#endif

   type :: crm_input_type

      real(crm_rknd), allocatable :: zmid(:,:)           ! Global grid height (m)
      real(crm_rknd), allocatable :: zint(:,:)           ! Global grid interface height (m)
      real(crm_rknd), allocatable :: tl(:,:)             ! Global grid temperature (K)
      real(crm_rknd), allocatable :: ql(:,:)             ! Global grid water vapor (g/g)
      real(crm_rknd), allocatable :: qccl(:,:)           ! Global grid cloud liquid water (g/g)
      real(crm_rknd), allocatable :: qiil(:,:)           ! Global grid cloud ice (g/g)
      real(crm_rknd), allocatable :: ps(:)               ! Global grid surface pressure (Pa)
      real(crm_rknd), allocatable :: pmid(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), allocatable :: pint(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), allocatable :: pdel(:,:)           ! Layer's pressure thickness (Pa)
      real(crm_rknd), allocatable :: phis(:)             ! Global grid surface geopotential (m2/s2)
      real(crm_rknd), allocatable :: ul(:,:)             ! Global grid u (m/s)
      real(crm_rknd), allocatable :: vl(:,:)             ! Global grid v (m/s)
      real(crm_rknd), allocatable :: ocnfrac(:)          ! area fraction of the ocean
      real(crm_rknd), allocatable :: tau00  (:)          ! large-scale surface stress (N/m2)
      real(crm_rknd), allocatable :: wndls  (:)          ! large-scale surface wind (m/s)
      real(crm_rknd), allocatable :: bflxls (:)          ! large-scale surface buoyancy flux (K m/s)
      real(crm_rknd), allocatable :: fluxu00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), allocatable :: fluxv00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), allocatable :: fluxt00(:)          ! surface sensible heat fluxes [K Kg/ (m2 s)]
      real(crm_rknd), allocatable :: fluxq00(:)          ! surface latent heat fluxes [ kg/(m2 s)]

      real(crm_rknd), allocatable :: naermod (:,:,:)     ! Aerosol number concentration [/m3]
      real(crm_rknd), allocatable :: vaerosol(:,:,:)     ! aerosol volume concentration [m3/m3]
      real(crm_rknd), allocatable :: hygro   (:,:,:)     ! hygroscopicity of aerosol mode 

      real(crm_rknd), allocatable :: ul_esmt(:,:)        ! input u for ESMT
      real(crm_rknd), allocatable :: vl_esmt(:,:)        ! input v for ESMT

      real(crm_rknd), allocatable :: t_vt(:,:)           ! CRM input of variance used for forcing tendency
      real(crm_rknd), allocatable :: q_vt(:,:)           ! CRM input of variance used for forcing tendency

      real(crm_rknd), allocatable :: relvar(:,:)          !  cloud liquid relative variance
      real(crm_rknd), allocatable :: nccn_prescribed(:,:) ! ccn prescribed concentration
      real(crm_rknd), allocatable :: ni_activated(:,:) 
      real(crm_rknd), allocatable :: npccn(:,:)
      real(crm_rknd), allocatable :: t_prev(:,:)
      real(crm_rknd), allocatable :: qv_prev(:,:)
      real(crm_rknd), allocatable :: ast(:,:)

      ! variable for shoc
      real(crm_rknd), allocatable :: sl(:,:)  
      real(crm_rknd), allocatable :: zm(:,:)
      real(crm_rknd), allocatable :: omega(:,:)
      
      real(crm_rknd), allocatable :: shf(:)      ! Sensible heat flux
      real(crm_rknd), allocatable :: cflx(:)       ! Latent heat flux
      real(crm_rknd), allocatable :: wsx(:)        ! Surface meridional momentum flux
      real(crm_rknd), allocatable :: wsy(:)        ! Surface zonal momentum flux  

      real(crm_rknd), allocatable :: tke_zt(:,:)  ! turbulent kinetic energy, interface
      real(crm_rknd), allocatable :: wthv(:,:) ! buoyancy flux
      real(crm_rknd), allocatable :: tkh(:,:)
      real(crm_rknd), allocatable :: tk(:,:)
      real(crm_rknd), allocatable :: alst(:,:)       ! liquid stratiform cloud fraction             [fraction]
      real(crm_rknd), allocatable :: qtracers(:,:,:) ! tracers
   end type crm_input_type
   !------------------------------------------------------------------------------------------------

contains
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_input_type
   subroutine crm_input_initialize(input, ncrms, nlev, MMF_microphysics_scheme)
      type(crm_input_type), intent(inout) :: input
      integer, intent(in) :: ncrms, nlev
      character(len=*), intent(in) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      if (.not. allocated(input%zmid))     allocate(input%zmid(ncrms,nlev))
      if (.not. allocated(input%zint))     allocate(input%zint(ncrms,nlev+1))
      if (.not. allocated(input%tl))       allocate(input%tl(ncrms,nlev))
      if (.not. allocated(input%ql))       allocate(input%ql(ncrms,nlev))
      if (.not. allocated(input%qccl))     allocate(input%qccl(ncrms,nlev))
      if (.not. allocated(input%qiil))     allocate(input%qiil(ncrms,nlev))
      if (.not. allocated(input%ps))       allocate(input%ps(ncrms))
      if (.not. allocated(input%pmid))     allocate(input%pmid(ncrms,nlev))
      if (.not. allocated(input%pint))     allocate(input%pint(ncrms,nlev+1))
      if (.not. allocated(input%pdel))     allocate(input%pdel(ncrms,nlev))
      if (.not. allocated(input%phis))     allocate(input%phis(ncrms))
      if (.not. allocated(input%ul))       allocate(input%ul(ncrms,nlev))
      if (.not. allocated(input%vl))       allocate(input%vl(ncrms,nlev))
      if (.not. allocated(input%ocnfrac))  allocate(input%ocnfrac(ncrms))
      if (.not. allocated(input%tau00))    allocate(input%tau00(ncrms))
      if (.not. allocated(input%wndls))    allocate(input%wndls(ncrms))
      if (.not. allocated(input%bflxls))   allocate(input%bflxls(ncrms))
      if (.not. allocated(input%fluxu00))  allocate(input%fluxu00(ncrms))
      if (.not. allocated(input%fluxv00))  allocate(input%fluxv00(ncrms))
      if (.not. allocated(input%fluxt00))  allocate(input%fluxt00(ncrms))
      if (.not. allocated(input%fluxq00))  allocate(input%fluxq00(ncrms))

      call prefetch(input%zmid)
      call prefetch(input%zint)
      call prefetch(input%tl)
      call prefetch(input%ql)
      call prefetch(input%qccl)
      call prefetch(input%qiil)
      call prefetch(input%ps)
      call prefetch(input%pmid)
      call prefetch(input%pint)
      call prefetch(input%pdel)
      call prefetch(input%phis)
      call prefetch(input%ul)
      call prefetch(input%vl)
      call prefetch(input%ocnfrac)
      call prefetch(input%tau00)
      call prefetch(input%wndls)
      call prefetch(input%bflxls)
      call prefetch(input%fluxu00)
      call prefetch(input%fluxv00)
      call prefetch(input%fluxt00)
      call prefetch(input%fluxq00)

      if (trim(MMF_microphysics_scheme) .eq. 'm2005') then
         if (.not. allocated(input%naermod))  allocate(input%naermod(ncrms,nlev,ntot_amode))
         if (.not. allocated(input%vaerosol)) allocate(input%vaerosol(ncrms,nlev,ntot_amode))
         if (.not. allocated(input%hygro))    allocate(input%hygro(ncrms,nlev,ntot_amode))
         call prefetch(input%naermod)
         call prefetch(input%vaerosol)
         call prefetch(input%hygro)
      end if

#if defined(MMF_ESMT)
      if (.not. allocated(input%ul_esmt))  allocate(input%ul_esmt(ncrms,nlev))
      if (.not. allocated(input%vl_esmt))  allocate(input%vl_esmt(ncrms,nlev))
#endif

      if (.not. allocated(input%t_vt)) allocate(input%t_vt(ncrms,nlev))
      if (.not. allocated(input%q_vt)) allocate(input%q_vt(ncrms,nlev))
      call prefetch(input%t_vt)
      call prefetch(input%q_vt)

      if (.not. allocated(input%relvar)) allocate(input%relvar(ncrms,nlev))
      if (.not. allocated(input%nccn_prescribed)) allocate(input%nccn_prescribed(ncrms, nlev))
      if (.not. allocated(input%npccn)) allocate(input%npccn(ncrms, nlev))
      if (.not. allocated(input%ni_activated)) allocate(input%ni_activated(ncrms, nlev))
      if (.not. allocated(input%t_prev)) allocate(input%t_prev(ncrms, nlev))
      if (.not. allocated(input%qv_prev)) allocate(input%qv_prev(ncrms, nlev))
      if (.not. allocated(input%sl)) allocate(input%sl(ncrms, nlev))
      if (.not. allocated(input%ast)) allocate(input%ast(ncrms, nlev))
      if (.not. allocated(input%omega)) allocate(input%omega(ncrms, nlev))
      if (.not. allocated(input%zm)) allocate(input%zm(ncrms, nlev))
      if (.not. allocated(input%shf)) allocate(input%shf(ncrms))
      if (.not. allocated(input%cflx)) allocate(input%cflx(ncrms))  
      if (.not. allocated(input%wsx)) allocate(input%wsx(ncrms))
      if (.not. allocated(input%wsy)) allocate(input%wsy(ncrms))

      if (.not. allocated(input%tke_zt)) allocate(input%tke_zt(ncrms,nlev))
      if (.not. allocated(input%wthv)) allocate(input%wthv(ncrms,nlev))
      if (.not. allocated(input%tkh)) allocate(input%tkh(ncrms,nlev))
      if (.not. allocated(input%tk)) allocate(input%tk(ncrms,nlev))
      if (.not. allocated(input%alst)) allocate(input%alst(ncrms,nlev))
      if (.not. allocated(input%qtracers)) allocate(input%qtracers(ncrms,nlev,10))

      call prefetch(input%relvar)
      call prefetch(input%nccn_prescribed)
      call prefetch(input%npccn)
      call prefetch(input%ni_activated)
      call prefetch(input%t_prev)
      call prefetch(input%qv_prev)
      call prefetch(input%sl)
      call prefetch(input%ast)
      call prefetch(input%omega)
      call prefetch(input%zm)
      call prefetch(input%shf)
      call prefetch(input%cflx)
      call prefetch(input%wsx)
      call prefetch(input%wsy)
      call prefetch(input%tke_zt)
      call prefetch(input%wthv)
      call prefetch(input%tkh)
      call prefetch(input%tk)
      call prefetch(input%alst)
      call prefetch(input%qtracers)

      ! Initialize
      input%zmid    = 0
      input%zint    = 0
      input%tl      = 0
      input%ql      = 0
      input%qccl    = 0
      input%qiil    = 0
      input%ps      = 0
      input%pmid    = 0
      input%pint    = 0
      input%pdel    = 0
      input%phis    = 0
      input%ul      = 0
      input%vl      = 0

      input%ocnfrac = 0
      input%tau00   = 0
      input%wndls   = 0
      input%bflxls  = 0
      input%fluxu00 = 0
      input%fluxv00 = 0
      input%fluxt00 = 0
      input%fluxq00 = 0

      input%relvar  = 0
      input%nccn_prescribed = 0
      input%npccn = 0
      input%t_prev = 0
      input%qv_prev = 0
      input%sl = 0
      input%ast = 0
      input%omega = 0
      input%zm = 0
      input%shf = 0
      input%cflx = 0
      input%wsx  = 0
      input%wsy  = 0
      input%tke_zt = 0
      input%wthv = 0
      input%tkh = 0
      input%tk = 0
      input%alst = 0

      if (trim(MMF_microphysics_scheme) .eq. 'm2005') then
         input%naermod  = 0
         input%vaerosol = 0
         input%hygro    = 0
      end if

#if defined( MMF_ESMT )
      input%ul_esmt = 0
      input%vl_esmt = 0
#endif

      input%t_vt = 0
      input%q_vt = 0

      input%relvar = 0
      input%nccn_prescribed = 0
   end subroutine crm_input_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_input_finalize(input, MMF_microphysics_scheme)
      type(crm_input_type), intent(inout) :: input
      character(len=*), intent(in) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      if (allocated(input%zmid))    deallocate(input%zmid)
      if (allocated(input%zint))    deallocate(input%zint)
      if (allocated(input%tl))      deallocate(input%tl)
      if (allocated(input%ql))      deallocate(input%ql)
      if (allocated(input%qccl))    deallocate(input%qccl)
      if (allocated(input%qiil))    deallocate(input%qiil)
      if (allocated(input%ps))      deallocate(input%ps)
      if (allocated(input%pmid))    deallocate(input%pmid)
      if (allocated(input%pint))    deallocate(input%pint)
      if (allocated(input%pdel))    deallocate(input%pdel)
      if (allocated(input%phis))    deallocate(input%phis)
      if (allocated(input%ul))      deallocate(input%ul)
      if (allocated(input%vl))      deallocate(input%vl)

      if (allocated(input%ocnfrac)) deallocate(input%ocnfrac)
      if (allocated(input%tau00))   deallocate(input%tau00)
      if (allocated(input%wndls))   deallocate(input%wndls)
      if (allocated(input%bflxls))  deallocate(input%bflxls)
      if (allocated(input%fluxu00)) deallocate(input%fluxu00)
      if (allocated(input%fluxv00)) deallocate(input%fluxv00)
      if (allocated(input%fluxt00)) deallocate(input%fluxt00)
      if (allocated(input%fluxq00)) deallocate(input%fluxq00)

      if (trim(MMF_microphysics_scheme) .eq. 'm2005') then
         if (allocated(input%naermod))    deallocate(input%naermod)
         if (allocated(input%vaerosol))   deallocate(input%vaerosol)
         if (allocated(input%hygro))      deallocate(input%hygro)
      end if

#if defined(MMF_ESMT)
      if (allocated(input%ul_esmt)) deallocate(input%ul_esmt)
      if (allocated(input%vl_esmt)) deallocate(input%vl_esmt)
#endif

      if (allocated(input%t_vt)) deallocate(input%t_vt)
      if (allocated(input%q_vt)) deallocate(input%q_vt)

      if (allocated(input%relvar)) deallocate(input%relvar)
      if (allocated(input%nccn_prescribed)) deallocate(input%nccn_prescribed)
      if (allocated(input%npccn)) deallocate(input%npccn)
      if (allocated(input%ni_activated)) deallocate(input%ni_activated)
      if (allocated(input%t_prev)) deallocate(input%t_prev)
      if (allocated(input%qv_prev)) deallocate(input%qv_prev)
      if (allocated(input%sl)) deallocate(input%sl)
      if (allocated(input%ast)) deallocate(input%ast)
      if (allocated(input%zm)) deallocate(input%zm)
      if (allocated(input%omega)) deallocate(input%omega)
      if (allocated(input%shf)) deallocate(input%shf)      
      if (allocated(input%cflx)) deallocate(input%cflx)       
      if (allocated(input%wsx)) deallocate(input%wsx)        
      if (allocated(input%wsy)) deallocate(input%wsy)

      if (allocated(input%tke_zt)) deallocate(input%tke_zt)
      if (allocated(input%wthv)) deallocate(input%wthv)
      if (allocated(input%tkh)) deallocate(input%tkh)
      if (allocated(input%tk)) deallocate(input%tk)
      if (allocated(input%alst)) deallocate(input%alst)
      if (allocated(input%qtracers)) deallocate(input%qtracers)
   end subroutine crm_input_finalize 

end module crm_input_module
