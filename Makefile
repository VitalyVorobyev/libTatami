Sources=cnvl.cpp conv_coef.cpp Rk.cpp RkRdetRnpPdf.cpp RbkgPdf.cpp ResConst.cpp Faddeeva.cpp parmanager.cpp wtag.cpp icpvevent.cpp toypdf.cpp toypdfgen.cpp abspdf.cpp absicpvpdf.cpp rrecpdf.cpp rascpdf.cpp ttools.cpp rascrnppars.cpp rrecpars.cpp rkparam.cpp ResVar.cpp
Executable=libtatami.so
CFlags=-c -Wall -fPIC -g -Iinc -std=c++14
LDFlags= -shared -fPIC -std=c++14 -I. -lm -lstdc++ 
ObjectDir=obj/
SourceDir=src/
BinDir=bin/

CC=g++
RM=rm

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
Objects=$(Sources:.cpp=.o)
CSources=$(addprefix $(SourceDir),$(Sources))
CObjects=$(addprefix $(ObjectDir),$(Objects))
CExecutable=$(addprefix $(BinDir),$(Executable))

all: $(CSources) $(CExecutable)

$(CExecutable): $(CObjects)
	$(CC) $(LDFlags) $(CObjects) -o $@

$(ObjectDir)%.o: $(SourceDir)%.cpp
	$(CC) $(CFlags) $< -o $@

clean:
	$(RM) $(CObjects)
