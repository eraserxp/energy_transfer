!==========================================================================================
      module molecule_information
! to save the values for dipole, electric_field, rot_const, lattice_constant
        implicit none
        double precision, save :: dipole
        double precision, save :: electric_field
        double precision, save :: rot_const
        double precision, save :: lattice_constant
        double precision, save :: alpha_parallel 
        double precision, save :: alpha_perpendicular
      end module molecule_information

      module some_parameters
        implicit none
        double precision, save :: time
        double precision, save :: phase 
        double precision, save :: t_var != 3*1.D-6
        integer, save :: n_iterate != 1000
        double precision, save :: time_interval
        integer, save :: imin
        integer, save :: imax
            
      end module some_parameters




!=======================================================================================
      subroutine initialize_molecule_information
        use molecule_information
        use some_parameters
        character(32) :: item
        ! for LiCs molecule
      	!rot_const=11.7998D9/2.D0  !unit=Hz 1GHz=10^9 Hz
      	!Dipole=5.529D0 
        !Lattice_Constant=4.D-7         
        !electric_field = 1.D5 ! V/m
        !alpha_parallel = 597.0D0 ! a.u.
        !alpha_perpendicular = 262.5D0 ! a.u.
        write(*,*) " "
        
        open(100,file='parameters_for_calculation.txt')
        read(100,*) item, rot_const
        write(*,*) item, " = ", rot_const
        
        read(100,*) item, Dipole
        write(*,*) item, " = ", Dipole
        
        read(100,*) item, Lattice_Constant
        write(*,*) item, " = ", Lattice_Constant
        
        read(100,*) item, electric_field
        write(*,*) item, " = ", electric_field
        
        read(100,*) item, alpha_parallel
        write(*,*) item, " = ", alpha_parallel
        
        read(100,*) item, alpha_perpendicular
        write(*,*) item, " = ", alpha_perpendicular
        
        read(100,*) item, t_var
        write(*,*) item, " = ", t_var
        
        read(100,*) item, n_iterate
        write(*,*) item, " = ", n_iterate  

!        read(100,*) item, imin
!        write(*,*) item, " = ", imin
        
!        read(100,*) item, imax
!        write(*,*) item, " = ", imax                    
        close(100)
        
        CALL Dipole_in_atomic_unit(Dipole)
        CALL Length_in_Bohr(Lattice_Constant)
        CALL Hz_to_atomic_unit(rot_const)
        CALL Field_in_atomic_unit(electric_field)
        
        time_interval = t_var/n_iterate 
        write(*,*) " "
        write(*,*) "Calculation begins"
      end subroutine initialize_molecule_information 


