info_TEXINFOS = cometsuite.texi
AM_MAKEINFOHTMLFLAGS = --css-include=$(top_srcdir)/doc/base.css
man_MANS=rundynamics.1

rundynamics.1: ../src/rundynamics
	$(HELP2MAN) -N ../src/rundynamics --output=$@
