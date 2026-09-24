add_subdirectory(libs)
add_subdirectory(tests)

if(BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
endif(BUILD_STATIC)

if (RPATH)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif (RPATH)

add_executable(nrscope main.cc)
add_executable(nrscan scan_main.cc)

set(SRSUE_SOURCES srsue_phy srsue_stack srsue_upper srsue_mac srsue_rrc srslog system)
set(SRSRAN_SOURCES srsran_common srsran_mac srsran_phy srsran_rf_utils srsran_radio srsran_gtpu srsran_rlc srsran_pdcp rrc_asn1 srslog support system )

set(SRSUE_SOURCES ${SRSUE_SOURCES} srsue_nr_stack srsue_rrc_nr srsue_mac_nr)
set(SRSRAN_SOURCES ${SRSRAN_SOURCES} rrc_nr_asn1 ngap_nr_asn1)


target_link_libraries(nrscope   srsue_phy
                                srsran_common
                                srsran_phy
                                srsran_radio
                                rrc_nr_asn1
                                nrscope_libs
                                ${SRSUE_SOURCES} 
                                ${SRSRAN_SOURCES}
                                ${CMAKE_THREAD_LIBS_INIT}
                                ${Boost_LIBRARIES}
                                ${YAML_CPP_LIBRARIES}
                                "liquid"
                                "liquid")

target_link_libraries(nrscan    srsue_phy
                                srsran_common
                                srsran_phy
                                srsran_radio
                                rrc_nr_asn1
                                nrscope_libs
                                ${SRSUE_SOURCES} 
                                ${SRSRAN_SOURCES}
                                ${CMAKE_THREAD_LIBS_INIT}
                                ${Boost_LIBRARIES}
                                ${YAML_CPP_LIBRARIES}
                                "liquid"
                                "liquid")

# Make the working config visible next to the binaries. nrscope and nrscan read
# "config.yaml" from the current working directory, so without this you have to
# copy the file in by hand -- and a stale copy silently overrides your edits.
#
# This links rather than copies, so that editing nrscope/config/config.yaml
# takes effect on the next run (no rebuild), and so that rebuilding can never
# overwrite the config: a symlink has no content to clobber.
#
# config.yaml is the default (the lab srsRAN gNB). Other configs are linked too
# and picked by passing them as the first argument, e.g.
#   sudo ./nrscope config_sunrise.yaml
set(NRSCOPE_CONFIG_DIR ${PROJECT_SOURCE_DIR}/nrscope/config)
set(NRSCOPE_CONFIGS config.yaml config_sunrise.yaml config_x410_UC200.yaml)

set(NRSCOPE_CONFIG_LINK_COMMANDS)
foreach(cfg ${NRSCOPE_CONFIGS})
  list(APPEND NRSCOPE_CONFIG_LINK_COMMANDS
       COMMAND ${CMAKE_COMMAND}
               -DSRC=${NRSCOPE_CONFIG_DIR}/${cfg}
               -DDST=${CMAKE_CURRENT_BINARY_DIR}/${cfg}
               -P ${PROJECT_SOURCE_DIR}/cmake/modules/LinkNrscopeConfig.cmake)
endforeach()

add_custom_target(nrscope_config ALL
  ${NRSCOPE_CONFIG_LINK_COMMANDS}
  COMMENT "Linking ${NRSCOPE_CONFIGS} from ${NRSCOPE_CONFIG_DIR}"
  VERBATIM)

add_dependencies(nrscope nrscope_config)
add_dependencies(nrscan  nrscope_config)

if(PUSH_TO_GOOGLE)
  target_link_libraries(nrscope ${PYTHON_LIBRARIES})
  target_link_libraries(nrscan ${PYTHON_LIBRARIES})
endif()

if (RPATH)
  set_target_properties(nrscope PROPERTIES INSTALL_RPATH ".")
  set_target_properties(nrscan PROPERTIES INSTALL_RPATH ".")
endif (RPATH)

install(TARGETS nrscope DESTINATION ${RUNTIME_DIR})
install(TARGETS nrscan DESTINATION ${RUNTIME_DIR})
                              