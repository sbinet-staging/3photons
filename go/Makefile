## simple makefile to build a stand-alone Go program
.PHONY: all build clean

all: build

build:
	@go build -o 3photons main.go rand.go event.go spinor.go p4.go

clean:
	@/bin/rm ./3photons

run: build
	@./3photons

