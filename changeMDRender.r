options(rstudio.markdownToHTML =
	function(inputFile, outputFile) {
		require(markdown)
		htmlOptions <- markdownHTMLOptions(defaults=TRUE)
		htmlOptions <- htmlOptions[htmlOptions != "hard_wrap"]
		markdownToHTML(inputFile, outputFile, options = htmlOptions)
	}
) 

