SourceDir=src/
Sources=$(notdir $(wildcard $(SourceDir)*.cpp))
Executable=libtatami.so
CFlags=-c -Wall -fPIC -g -Iinc -I. `root-config --cflags` -std=c++14
LDFlags= -shared -fPIC -I. -lm -lstdc++ -lfaddeeva -Wl,--no-as-needed `root-config --glibs` -std=c++14
ObjectDir=obj/
IncludeDir=include/
BinDir=bin/

CC=g++ -O2
RM=rm

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
Objects=$(Sources:.cpp=.o)
CSources=$(addprefix $(SourceDir),$(Sources))
CObjects=$(addprefix $(ObjectDir),$(Objects))
CExecutable=$(addprefix $(BinDir),$(Executable))
CIncs = -I$(IncludeDir)

all: $(CSources) $(CExecutable)

$(CExecutable): $(CObjects)
	$(CC) $(LDFlags) $(CObjects) -o $@

$(ObjectDir)%.o: $(SourceDir)%.cpp
	$(CC) $(CFlags) $(CIncs) $< -o $@

clean:
	$(RM) $(CObjects)
