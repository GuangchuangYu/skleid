PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: docs check clean

docs:
	Rscript -e 'roxygen2::roxygenise(".")'

build:
	cd ..;\
	R CMD build $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

windows:
	Rscript -e 'rhub::check_on_windows(".")';\
	sleep 10;

windowsstatus:
	STATUS := $(shell Rscript -e 'ypages:::check_rhub_status()')

addtorepo: windows
	cd ../drat;\
	Rscript -e 'drat:::insert("../$(PKGNAME)_$(PKGVERS).tar.gz", "docs")';\
	Rscript -e 'drat:::insert(ypages::get_windows_binary(), "docs")';\
	git add .; git commit -m '$(PKGNAME)_$(PKGVERS)'; git push -u origin master
