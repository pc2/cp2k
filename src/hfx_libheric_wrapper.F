MODULE libheric_wrapper
#if defined(__LIBHERIC)
    USE libheric_interface, ONLY: hfx_set_quartet_class => accelerator_context_set_quartet_class
#endif
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: hfx_set_quartet_class

CONTAINS
#if !defined (__LIBHERIC)
    SUBROUTINE hfx_set_quartet_class(ctx, la, lb, lc, ld)
       INTEGER :: ctx, la, lb, lc, ld
    END SUBROUTINE
#endif
END MODULE libheric_wrapper
