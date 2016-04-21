# mk_diags

modele_set_dependencies()
mkdiags_set_dependencies()

# Sources
set(srcs_f
   agcstat.f
   csregrid.f
   diffreport.f
   miscnc.f
   prtadiurn.f
   prtaj.f
   prtajl.f
   prtareg.f
   prtconsrv.f
   prtisccp.f
   prtolnst.f
   prtostat.f
   prtotj.f
   prtrvr.f
   prtspeca.f
   scaleaccdrv.f
   scaleacc.f
   sumfiles.f
   write_2d_as_giss4d.f
   write_giss2d.f
)

#
# Loop over diagnostics categories and compile each Fortran print
# routine with a driver program which opens/closes the input file
# and passes any extra command-line arguments to the print routine
#
set(diags
   aj areg ajl adiurn consrv rvr isccp olnst otj
)
foreach(diag ${diags})
   set(PRTDIAG "prt${diag}")
   configure_file (
     "${CMAKE_CURRENT_SOURCE_DIR}/prtdrv.f.in"
     "${CMAKE_CURRENT_BINARY_DIR}/${PRTDIAG}drv.f"
   )
   add_executable(${PRTDIAG} 
      ${CMAKE_CURRENT_BINARY_DIR}/${PRTDIAG}drv.f 
      ${PRTDIAG}.f
      miscnc.f)
   #message(STATUS " ---> Add executable ${PRTDIAG}")
   target_link_libraries(${PRTDIAG}  ${MKDIAGS_EXTERNAL_LIBS})
   install(TARGETS ${PRTDIAG} DESTINATION bin)
endforeach()

#
# Compile the generic scaling routine with an appropriate driver
#
add_executable(scaleacc scaleaccdrv.f scaleacc.f csregrid.f miscnc.f)
target_link_libraries(scaleacc ${MKDIAGS_EXTERNAL_LIBS})
install(TARGETS scaleacc DESTINATION bin)

#
# Compile miscellaneous routines
#
set(apps
  agcstat sumfiles write_giss2d write_2d_as_giss4d diffreport prtostat prtspeca
)
foreach(app ${apps})
   add_executable(${app} ${app}.f miscnc.f)
   #message(STATUS " ---> Add executable ${app}")
   target_link_libraries(${app} ${ModelE_EXTERNAL_LIBS})
   install(TARGETS ${app} DESTINATION bin)
endforeach()

