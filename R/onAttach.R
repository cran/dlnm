.onAttach <- function(lib, pkg) {
	meta <- packageDescription("dlnm")
	attachmsg <- paste("This is dlnm ",meta$Version,
		". For details type: vignette('dlnmOverview').",sep="")
	packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}