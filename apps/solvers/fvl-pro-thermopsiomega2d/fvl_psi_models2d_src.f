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

subroutine psi_initialize_pseudo_patches()
      integer(kind=__fvl_integer_kind__)::i,numedges
      integer(kind=__fvl_integer_kind__),dimension(1:mesh%getnumedges())::edgeindices
      type(fvl_patch2d),pointer::patch
      if(psi_numpseudoboundconds>0) then
            numedges=0
            do i=1,psi_numpseudoboundconds
                  patch=>psi_pseudoboundconds(i)%getpatch()
                  edgeindices(numedges+1:numedges+patch%getnumedges())=patch%getedgeindices()
                  numedges=numedges+patch%getnumedges()
            end do
            psi_pseudo_patch=fvl_patch2d(mesh=mesh,vertindices=[integer::],edgeindices=edgeindices(1:numedges),cellindices=[integer::])
      else
            psi_pseudo_patch=fvl_patch2d(mesh=mesh,vertindices=[integer::],edgeindices=[integer::],cellindices=[integer::])
      end if
end subroutine psi_initialize_pseudo_patches

! ----------------------------------------------------------------------------
! inner fields
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! boundary fields
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! boundary conditions
! ----------------------------------------------------------------------------

subroutine fvl_set_psi_constpseudovalue_boundcond1(codes)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_constpseudovalue_boundcond2d_init(patch=psi_pseudopatches(psi_numpseudopatches),&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_constpseudovalue_boundcond1

subroutine fvl_set_psi_constpseudomixed_boundcond1(codes,value,coefs)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_constpseudomixed_boundcond2d_init1(patch=psi_pseudopatches(psi_numpseudopatches),value=value,coefs=coefs,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_constpseudomixed_boundcond1

subroutine fvl_set_psi_constpseudomixed_boundcond2(codes,value,coefsfun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_psi_coefsfun)::coefsfun
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_constpseudomixed_boundcond2d_init2(patch=psi_pseudopatches(psi_numpseudopatches),value=value,coefsfun=coefsfun,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_constpseudomixed_boundcond2

subroutine fvl_set_psi_fixedpseudomixed_boundcond1(codes,valuefun,coefs)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      procedure(fvl_psi_valuefun)::valuefun
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_fixedpseudomixed_boundcond2d_init1(patch=psi_pseudopatches(psi_numpseudopatches),valuefun=valuefun,coefs=coefs,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_fixedpseudomixed_boundcond1

subroutine fvl_set_psi_fixedpseudomixed_boundcond2(codes,value,coefsfun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_psi_coefsfun)::coefsfun
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_fixedpseudomixed_boundcond2d_init2(patch=psi_pseudopatches(psi_numpseudopatches),value=value,coefsfun=coefsfun,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_fixedpseudomixed_boundcond2

subroutine fvl_set_psi_fixedpseudomixed_boundcond3(codes,valuefun,coefsfun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      procedure(fvl_psi_valuefun)::valuefun
      procedure(fvl_psi_coefsfun)::coefsfun
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_fixedpseudomixed_boundcond2d_init3(patch=psi_pseudopatches(psi_numpseudopatches),valuefun=valuefun,coefsfun=coefsfun,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_fixedpseudomixed_boundcond3

subroutine fvl_set_psi_constpseudocombined_boundcond1(codes,value,coefs)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_constpseudocombined_boundcond2d_init1(patch=psi_pseudopatches(psi_numpseudopatches),value=value,coefs=coefs,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_constpseudocombined_boundcond1

subroutine fvl_set_psi_constpseudocombined_boundcond2(codes,value,coefsfun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_psi_coefsfun)::coefsfun
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_constpseudocombined_boundcond2d_init2(patch=psi_pseudopatches(psi_numpseudopatches),value=value,coefsfun=coefsfun,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_constpseudocombined_boundcond2

subroutine fvl_set_psi_fixedpseudocombined_boundcond1(codes,valuefun,coefs)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      procedure(fvl_psi_valuefun)::valuefun
      real(kind=__fvl_real_kind__),dimension(1:2),intent(in)::coefs
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_fixedpseudocombined_boundcond2d_init1(patch=psi_pseudopatches(psi_numpseudopatches),valuefun=valuefun,coefs=coefs,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_fixedpseudocombined_boundcond1

subroutine fvl_set_psi_fixedpseudocombined_boundcond2(codes,value,coefsfun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      real(kind=__fvl_real_kind__),intent(in)::value
      procedure(fvl_psi_coefsfun)::coefsfun
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_fixedpseudocombined_boundcond2d_init2(patch=psi_pseudopatches(psi_numpseudopatches),value=value,coefsfun=coefsfun,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_fixedpseudocombined_boundcond2

subroutine fvl_set_psi_fixedpseudocombined_boundcond3(codes,valuefun,coefsfun)
      integer(kind=__fvl_integer_kind__),dimension(1:),contiguous,intent(in)::codes
      procedure(fvl_psi_valuefun)::valuefun
      procedure(fvl_psi_coefsfun)::coefsfun
      ! create mesh patch
      psi_numpseudopatches=psi_numpseudopatches+1
      psi_pseudopatches(psi_numpseudopatches)=fvl_patch2d(mesh=mesh,codes=codes,edgespatch=.true.)
      ! crate single field
      psi_numsinglefields=psi_numsinglefields+1
      psi_singlefields(psi_numsinglefields)=fvl_scalarsinglefield2d(patch=psi_pseudopatches(psi_numpseudopatches),value=0.0d0,&
            label="psi"//fvl_trim(fvl_char(psi_numsinglefields)),filepath=solutionfilesdir,fileform=solutionfilesform)
      ! create boundary condition
      psi_numpseudoboundconds=psi_numpseudoboundconds+1
      psi_pseudoboundconds(psi_numpseudoboundconds)=fvl_psi_fixedpseudocombined_boundcond2d_init3(patch=psi_pseudopatches(psi_numpseudopatches),valuefun=valuefun,coefsfun=coefsfun,&
            psisinglefield=psi_singlefields(psi_numsinglefields),omegaconvinnerfield=omegaconv_innerfield,omegadiffinnerfield=omegadiff_innerfield,&
            forcemodel=fup_model)
      ! set boundary condition
      call fvl_set_psi_boundcond(psi_pseudoboundconds(psi_numpseudoboundconds))
end subroutine fvl_set_psi_fixedpseudocombined_boundcond3

! end of file
