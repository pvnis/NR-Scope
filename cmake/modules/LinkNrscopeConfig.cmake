# Links the working config.yaml into a target's runtime directory.
#
# nrscope and nrscan both open "config.yaml" relative to the current working
# directory, so a copy has to exist next to the binaries. A symlink is used
# rather than a copy on purpose:
#
#   * the link has no content of its own, so building can never overwrite or
#     wipe out the real config;
#   * edits to nrscope/config/config.yaml take effect on the next run, with
#     no rebuild needed.
#
# Invoked at build time via `cmake -P` with -DSRC=<source yaml> -DDST=<link>.

if(EXISTS "${DST}" AND NOT IS_SYMLINK "${DST}")
  # Somebody put a real file here (most likely an older manual copy). Never
  # delete it, but say so loudly: the binaries will read THAT file and silently
  # ignore nrscope/config/config.yaml.
  message(WARNING
    "${DST} is a regular file, not a link, so it is being left alone.\n"
    "   nrscope/nrscan will read it INSTEAD of ${SRC}.\n"
    "   Delete it to have the config linked automatically.")
else()
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRC}" "${DST}")
endif()
