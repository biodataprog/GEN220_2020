MARKDOWN=pandoc --from markdown --standalone --indented-code-classes=python

all: pdfall

pdfall: $(patsubst %.md,%.pdf,$(wildcard *.md)) 

htmlall: $(patsubst %.md,%.html,$(wildcard *.md)) 

%.html: %.md
	$(MARKDOWN) $< --to html5 --output $@ metadata.yaml

%.pdf: %.md
	$(MARKDOWN) $< --to latex --output $@ metadata.yaml
