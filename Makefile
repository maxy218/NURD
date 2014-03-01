SRCS = main.cpp algorithm.cpp class.cpp common.cpp const.cpp 
OBJS = $(SRCS:.c=.cpp)
LIBS = ./boost/boost_1_55
EXECUTABLE = NURD

all:
	g++ -O3 -o $(EXECUTABLE) -I $(LIBS) $(SRCS) 

.PHONY: clean

clean:
	@echo "cleaning project"
	-rm $(OBJ_DIR)/*.o
	@echo "clean completed"
