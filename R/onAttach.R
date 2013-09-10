###
### R routines for the R package dlnm (c) Antonio Gasparrini 2012-2013
#
`.onAttach` <- 
function(lib, pkg) {
#
################################################################################
#
  meta <- packageDescription("dlnm")
  attachmsg <- paste("This is dlnm ",meta$Version,
    ". For details: help(dlnm) and vignette('dlnmOverview').","\n",
    "Important changes since version 2.0.0","\n",
    "See: 'file.show(system.file('Changesince200',package='dlnm'))'",sep="")
  packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}

