!******************Convert Debye to atomic units*************************
SUBROUTINE Dipole_in_atomic_unit(Debye)
  IMPLICIT NONE
  DOUBLE PRECISION :: Debye
  Debye=Debye*0.3934302014076827
END SUBROUTINE Dipole_in_atomic_unit

function debye_to_au(Debye)
  IMPLICIT NONE
  double precision :: debye_to_au
  DOUBLE PRECISION :: Debye
  debye_to_au=Debye*0.3934302014076827
  return
END function debye_to_au




!*****************Convert meter to atomic unit****************************
SUBROUTINE Length_in_Bohr(meter)
  IMPLICIT NONE
  DOUBLE PRECISION :: meter
  meter=meter*1.889726133921252D10
END SUBROUTINE Length_in_Bohr	

function meter_to_au(meter)
  IMPLICIT NONE
  double precision :: meter_to_au
  DOUBLE PRECISION :: meter
  meter_to_au=meter*1.889726133921252D10
  return
END function meter_to_au




!*************Convert from energy atomic unit  to kHz**********************
SUBROUTINE Energy_in_kHz(energy_in_atomic_unit)
  IMPLICIT NONE
  DOUBLE PRECISION :: energy_in_atomic_unit
  energy_in_atomic_unit=energy_in_atomic_unit*6.57968392072144D12
END SUBROUTINE Energy_in_kHz    

function hartree_to_kHz(energy_in_atomic_unit)
  IMPLICIT NONE
  DOUBLE PRECISION :: energy_in_atomic_unit
  double precision :: hartree_to_kHz
  hartree_to_kHz=energy_in_atomic_unit*6.57968392072144D12
  return
END function hartree_to_kHz




!*************Convert from energy atomic unit  to MHz**********************
SUBROUTINE Energy_in_MHz(energy_in_atomic_unit)
  IMPLICIT NONE
  DOUBLE PRECISION :: energy_in_atomic_unit
  energy_in_atomic_unit=energy_in_atomic_unit*6.57968392072144D9
END SUBROUTINE Energy_in_MHz
 
 
 
 
!*****************Convert electric field to atomic unit****************************
SUBROUTINE Field_in_atomic_unit(electric_field)
  IMPLICIT NONE
  DOUBLE PRECISION :: electric_field !V/m
  electric_field=electric_field*1.944690567144141D-12
END SUBROUTINE Field_in_atomic_unit
 
function V_m_to_au(electric_field)
  IMPLICIT NONE
  double precision :: V_m_to_au
  DOUBLE PRECISION :: electric_field
  V_m_to_au=electric_field*1.944690567144141D-12
END function V_m_to_au
 






!*****************Convert Hz to atomic unit****************************
SUBROUTINE Hz_to_atomic_unit(hz)
  IMPLICIT NONE
  DOUBLE PRECISION :: hz
  hz=hz*2.418884324306202D-17
END SUBROUTINE Hz_to_atomic_unit	

function Hz_to_au(hz)
  IMPLICIT NONE
  double precision :: Hz_to_au
  DOUBLE PRECISION :: hz
  Hz_to_au=Hz*2.418884324306202D-17
  return
END function Hz_to_au




!*****************Convert second to atomic unit*************************
SUBROUTINE time_in_atomic_unit(time)
  IMPLICIT NONE
  DOUBLE PRECISION :: time
  time = time*4.1341373374D16
END SUBROUTINE time_in_atomic_unit	

function second_to_au(time)
  IMPLICIT NONE
  double precision :: second_to_au
  DOUBLE PRECISION :: time
  second_to_au = time*4.1341373374D16
  return
END function second_to_au



!****************Convert intensity to atomic unit***********************
subroutine intensity_in_au(intensity)
  implicit none
  double precision :: intensity ! in W/cm^2  
  intensity = intensity/3.50944758D16
end subroutine intensity_in_au


Function w_cm2_to_au(intensity)
  implicit none
  double precision :: w_cm2_to_au
  double precision :: intensity ! in W/cm^2  
  w_cm2_to_au = intensity/3.50944758D16
end function w_cm2_to_au
