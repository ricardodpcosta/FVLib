! _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!|  _ _ _   _   _   _      |                                   |
!| |_|_|_| |_| |_| |_|     |  The Great Finite Volume Library  |
!| |_|_    |_|_|_| |_|     |  Author: Ricardo Costa            |
!| |_|_|    \_\_/  |_|_ _  |  Version: 1.0                     |
!| |_|       \_/   |_|_|_| |  Release: January, 2022           |
!|_ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|
! About: three-dimensional fluid models
! Modification: April, 2024

! ============================================================================
! FUNCTIONS
! ============================================================================

! ----------------------------------------------------------------------------
! patches
! ----------------------------------------------------------------------------

subroutine omega_initialize_pseudo_patches()
      integer(kind=__fvl_integer_kind__)::i,numedges
      integer(kind=__fvl_integer_kind__),dimension(1:mesh%getnumedges())::edgeindices
      type(fvl_patch2d),pointer::patch
      if(omega_numpseudoboundconds>0) then
            numedges=0
            do i=1,omega_numpseudoboundconds
                  patch=>omega_pseudoboundconds(i)%getpatch()
                  edgeindices(numedges+1:numedges+patch%getnumedges())=patch%getedgeindices()
                  numedges=numedges+patch%getnumedges()
            end do
            omega_pseudo_patch=fvl_patch2d(mesh=mesh,vertindices=[integer::],edgeindices=edgeindices(1:numedges),cellindices=[integer::])
      else
            omega_pseudo_patch=fvl_patch2d(mesh=mesh,vertindices=[integer::],edgeindices=[integer::],cellindices=[integer::])
      end if
end subroutine omega_initialize_pseudo_patches

! ----------------------------------------------------------------------------
! inner fields
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! boundary fields
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! boundary conditions
! ----------------------------------------------------------------------------

subroutine fvl_set_omega_constpseudovalue_boundcond(codes,value)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      ! create mesh patch
      omega_numpseudopatches=omega_numpseudopatches+1
      omega_pseudopatches(omega_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! create boundary condition
      omega_numpseudoboundconds=omega_numpseudoboundconds+1
      omega_pseudoboundconds(omega_numpseudoboundconds)=fvl_omega_constpseudovalue_boundcond2d_init(patch=omega_pseudopatches(omega_numpseudopatches),value=value,&
            psinnboundfield=psinn_boundfield)
      ! set boundary condition
      call fvl_set_omega_boundcond(omega_pseudoboundconds(omega_numpseudoboundconds))
end subroutine fvl_set_omega_constpseudovalue_boundcond

subroutine fvl_set_omega_fixedpseudovalue_boundcond(codes,valuefun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      procedure(fvl_omega_valuefun)::valuefun
      ! create mesh patch
      omega_numpseudopatches=omega_numpseudopatches+1
      omega_pseudopatches(omega_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! create boundary condition
      omega_numpseudoboundconds=omega_numpseudoboundconds+1
      omega_pseudoboundconds(omega_numpseudoboundconds)=fvl_omega_fixedpseudovalue_boundcond2d_init(patch=omega_pseudopatches(omega_numpseudopatches),valuefun=valuefun,&
            psinnboundfield=psinn_boundfield)
      ! set boundary condition
      call fvl_set_omega_boundcond(omega_pseudoboundconds(omega_numpseudoboundconds))
end subroutine fvl_set_omega_fixedpseudovalue_boundcond

! end of file
