# generate a latex compatible version of the html for use by pandoc and pdflatex

require(markdown)
htmlOptions <- markdownHTMLOptions(defaults=TRUE)
htmlOptions <- htmlOptions[!(htmlOptions %in% c("hard_wrap", "base64_images")) ]
markdownToHTML("flightetal_draft.md", "flightetal_draft_noImages.html", options = htmlOptions, stylesheet='htmlpub.css')

system("pandoc flightetal_draft_noImages.html -o flightetal_draft.pdf")