      subroutine vtune_resume
        external ftn_vtune_resume
        write (*,*) 'VTune: resuming collection'
        call ftn_vtune_resume()
      end subroutine

      subroutine vtune_pause
        external ftn_vtune_pause
        write (*,*) 'VTune: pausing collection'
        call ftn_vtune_pause()
      end subroutine

