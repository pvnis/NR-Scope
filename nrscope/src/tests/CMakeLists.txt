file(GLOB_RECURSE SOURCES "*.cc" "*.c")
add_library(test_libs STATIC ${SOURCES})

# message(${PROJECT_SOURCE_DIR})
# target_include_directories(ngscope_dci PUBLIC ${PROJECT_SOURCE_DIR}/ngscope/hdr)

if(ENABLE_GUI AND SRSGUI_FOUND)
  target_link_libraries(test_libs ${SRSGUI_LIBRARIES})

target_link_libraries(test_libs ${SRSRAN_SOURCES}
                                   ${SRSUE_SOURCES})
endif()

add_executable(test_comm test_comm_path.cc)
add_executable(google_client google_client_test.cc)
add_executable(google_client_cc google_client_cc_test.cc)
add_executable(rrc_reconf_test rrc_reconf_test.cc)

# Ground-truth test of the sensing pipeline: synthesises a channel from known
# paths and checks the range-Doppler map puts them where they were placed. Needs
# no radio and no gNB, so it can run on any machine and in CI.
add_executable(sensing_selftest sensing_selftest.c)
target_link_libraries(sensing_selftest srsue_phy
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
                                       m)


target_link_libraries(test_comm srsue_phy
                                srsran_common
                                srsran_phy
                                srsran_radio
                                rrc_nr_asn1
                                nrscope_libs
                                test_libs
                                ${SRSUE_SOURCES} 
                                ${SRSRAN_SOURCES}
                                ${CMAKE_THREAD_LIBS_INIT}
                                ${Boost_LIBRARIES}
                                ${YAML_CPP_LIBRARIES})

target_link_libraries(google_client_cc srsue_phy
                                       srsran_common
                                       srsran_phy
                                       srsran_radio
                                       rrc_nr_asn1
                                       nrscope_libs
                                       test_libs
                                       ${SRSUE_SOURCES} 
                                       ${SRSRAN_SOURCES}
                                       ${CMAKE_THREAD_LIBS_INIT}
                                       ${Boost_LIBRARIES}
                                       ${YAML_CPP_LIBRARIES})    
                                       
target_link_libraries(rrc_reconf_test srsue_phy
                                      srsran_common
                                      srsran_phy
                                      srsran_radio
                                      rrc_nr_asn1
                                      nrscope_libs
                                      test_libs
                                      ${SRSUE_SOURCES} 
                                      ${SRSRAN_SOURCES}
                                      ${CMAKE_THREAD_LIBS_INIT}
                                      ${Boost_LIBRARIES}
                                      ${YAML_CPP_LIBRARIES}) 

if(PUSH_TO_GOOGLE)
  target_link_libraries(google_client ${PYTHON_LIBRARIES})
  # target_link_libraries(google_client_cc ${PYTHON_LIBRARIES})
  # target_link_libraries(google_client_cc ${GOOGLE_CLOUD_CPP_STRORAGE_LIBRARIES})
endif()


# NR-Scope's own PDSCH DM-RS regeneration against srsRAN's mapper; see the file.
# Sensing divides received DM-RS by these pilots, so this must keep passing.
add_executable(dmrs_check_test dmrs_check_test.cc)
target_link_libraries(dmrs_check_test nrscope_libs
                                      srsran_phy
                                      srsran_common
                                      ${CMAKE_THREAD_LIBS_INIT})
add_test(nrscope_dmrs_check dmrs_check_test)

# Target detection under the TDD sampling pattern (nr_ue_tdd_detect): the replica
# the pattern paints is stronger than the target on the map, and must not win.
add_executable(tdd_detect_test tdd_detect_test.c)
target_link_libraries(tdd_detect_test nrscope_libs
                                      srsran_phy
                                      srsran_common
                                      ${CMAKE_THREAD_LIBS_INIT}
                                      m)
add_test(nrscope_tdd_detect tdd_detect_test)
add_test(nrscope_sensing_selftest sensing_selftest)

# Replays a sensing.record_estimates file through the sensing pipeline, offline: the
# recorded run's channel estimates, impairments included, through the live code.
add_executable(sensing_offline sensing_offline.cc)
target_link_libraries(sensing_offline srsue_phy
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
                                      m)
